import os
import time
import itertools
from functools import reduce

import sqlite3
from natsort import natsorted

from tables_functions import create_links, get_data, get_table_info

def isTableVisible(tableInfo):
    return tableInfo['options']['visible'] 

def getConstraintTablesConstraintKeyName(tablesInfo):
    return list(map(lambda t: {'tableName': t['tableName'], 'constraintKeyName': t['options']['pk']}, filter(isTableVisible, tablesInfo)))

class Linker:
    def __init__(self, tablesInfo, dataset_id, path, do_init=False):
        self.tablesInfo = tablesInfo

        # use this to get consistent ordering of columns, can't rely on dicts
        self.tableColumns = {t['tableName']: [x for x in t['tableData'][0].keys()] for t in tablesInfo}
        self.columnIndices = {t['tableName']: {x: i for i, x in enumerate(t['tableData'][0].keys())} for t in tablesInfo}
        self.tablesMap = {t['tableName']: t for t in tablesInfo}

        self.defaultConstraints = self.getDefaultConstraints()
        self.fieldNames = self.getFieldNames()
        self.tableIdToIdColumnMap = self.getTableIdToIdColumnMap()

        self.selections = self.emptySelections()
        self.numSelected = self.countNumSelected()
        self.totalSelected = self.getTotalSelected()
        self.queryResult = None
        self.sqlManager = SqlManager(self.tablesInfo, path, do_init)

    def getDefaultConstraints(self):
        ctables = getConstraintTablesConstraintKeyName(self.tablesInfo)
        # this is just a list of dicts with {'tableName': ..., 'constraintKeyName': ...} entries 
        # giving the primary keys of each table
        constraints = {}
        for ct in ctables:
            constraints[ct['tableName']] = self.getKeys(self.tablesInfo, ct['tableName'], ct['constraintKeyName'])
 
        return constraints

    def getFieldNames(self):
        return list(map(lambda ti: {'tableName': ti['tableName'], 'fieldNames': list(ti['tableData'][0].keys())}, filter(isTableVisible, self.tablesInfo)))

    def getTableIdToIdColumnMap(self):
        return {t['tableName']: t['constraintKeyName'] for t in getConstraintTablesConstraintKeyName(self.tablesInfo)}

    def getKeys(self, tablesInfo, tableName, k):
        # get the values of the key column of each row for the given table
        table = self.tablesMap[tableName]
        if not isTableVisible(table):
            return []

        return [row[k] for row in table['tableData']]

    def countNumSelected(self):
        temp = getConstraintTablesConstraintKeyName(self.tablesInfo)
        results = {}
        for tableInfo in temp:
            tname = tableInfo['tableName']
            results[tname] = len(self.selections[tname])
        return results

    def getTotalSelected(self):
        return reduce(lambda a, b: a + b, self.numSelected.values())

    def emptySelections(self):
        temp = getConstraintTablesConstraintKeyName(self.tablesInfo)
        results = {}
        for tableInfo in temp:
            tname = tableInfo['tableName']
            results[tname] = []
            
        return results

    def getTableNames(self, tablesInfo):
        # (my addition) return names of visible tables
        visible_tables = filter(isTableVisible, tablesInfo)
        return [t['tableName'] for t in visible_tables]

    def query(self):
        fieldNames = self.fieldNames
        tableNames = self.getTableNames(self.tablesInfo)
        tableData = {tableName: [] for tableName in tableNames}

        # special case here if nothing is selected
        if self.getTotalSelected() == 0:
            # TODO could also check self.numSelected[tableName] == 0 to do individually

            # normally this joins all the tables together, runs a query and then 
            # splits the results back up into the individual tables for display. 
            # however when you just want to clear/reset them to their original state
            # it does a ton of useless work as the data required is already there in
            # exactly the right format 
            print('doing optimization of query')
            for tableName in tableNames:
                tableData[tableName] = self.sqlManager.simpleQuery(tableName, self.tableColumns[tableName])

            self.queryResult = tableData
            return

        constraints = self.makeConstraints()

        # this is a bit different from the original, because unlike alasql you
        # can't perform a query and then run other queries on the result of that
        # query (the equivalent in the Python/SQLite API would be running a query
        # on a cursor). So instead what this does is:
        # - construct and execute the same initial query as the original code uses...
        # - but wrap it into a "CREATE TEMP TABLES AS (<query>)" statement
        # - then run the other queries against this temporary table
        t = time.time()
        self.sqlManager.createTempTable(self.tablesInfo, constraints)
        print('createTempTable (totSel={}) took {:.2f}s'.format(self.getTotalSelected(), time.time() - t))

        # now do the equivalent of the prefixQuery method from the original, by 
        # running a SELECT DISTINCT query against the temp table, once per visible table

        for tableName in tableNames:
            tableData[tableName] = self.sqlManager.prefixQuery(tableName, self.tableColumns[tableName])

        self.queryResult = tableData

        # dispose of the temp table
        self.sqlManager.deleteTempTable()

    def addConstraint(self, tableName, rowIndex, dataSource):
        idColumn = self.tableIdToIdColumnMap[tableName]
        idVal = dataSource.data[idColumn][rowIndex]
        self.selections[tableName].append({'idVal': idVal, rowIndex: rowIndex})
        self.numSelected = self.countNumSelected()
        self.totalSelected = self.getTotalSelected()

    def removeConstraints(self, tableName):
        self.selections[tableName] = []
        self.numSelected = self.countNumSelected()
        self.totalSelected = self.getTotalSelected()

    def makeConstraints(self):
        temp = getConstraintTablesConstraintKeyName(self.tablesInfo)
        results = {}
        for tableInfo in temp:
            tname = tableInfo['tableName']
            results[tname] = self.selectionToConstraint(tname)
        return results

    def selectionToConstraint(self, tableName):
        if self.numSelected[tableName] == 0:
            return self.defaultConstraints[tableName]
        else:
            return list(map(lambda x: x['idVal'], self.selections[tableName]))

    def updateDataSources(self, dataSources):
        if self.totalSelected > 0:
            fieldNames = self.fieldNames
            queryResult = self.queryResult

            for i in range(len(fieldNames)):
                tableFieldName = fieldNames[i]
                tableName = tableFieldName['tableName']
                if self.numSelected[tableName] == 0:
                    dataSource = dataSources[tableName]
                    self.updateSingleTable(tableName, queryResult, dataSource)
        else:
            self.resetTables(dataSources)

    def updateSingleTable(self, tableName, queryResult, dataSource):
        newData = {}
        tableQueryResult = queryResult[tableName]

        idColumn = self.tableIdToIdColumnMap[tableName]
        # natural sorting on alphanumerical strings
        queryResult = natsorted(tableQueryResult, key=lambda x: x[idColumn])

        # for i, row in enumerate(queryResult[tableName]):
        for i, row in enumerate(queryResult):
            for prop in row:
                value = row[prop]
                if prop not in newData:
                    newData[prop] = []
                newData[prop].append(value)

        dataSource.data = newData

    def resetTables(self, dataSources):
        fieldNames = self.fieldNames
        queryResult = self.queryResult

        for i in range(len(fieldNames)):
            tableFieldName = fieldNames[i]
            tableName = tableFieldName['tableName']
            dataSource = dataSources[tableName]
            self.updateSingleTable(tableName, queryResult, dataSource)


class SqlManager:

    def __init__(self, tablesInfo, path='linker.sqlite3', do_init=False):
        self.db = sqlite3.connect(path)
        self.db.row_factory = sqlite3.Row

        # useful for debugging
        # self.db.set_trace_callback(print)

        # use this to get consistent ordering of columns, can't rely on dicts
        self.tableColumns = {t['tableName']: [x for x in t['tableData'][0].keys()] for t in tablesInfo}
        self.columnIndices = {t['tableName']: {x: i for i, x in enumerate(t['tableData'][0].keys())} for t in tablesInfo}

        if do_init:
            print('SqlManager: Doing tables init')
            self.initialiseTables(tablesInfo)
        self.firstTable = self.getFirstTable(tablesInfo)
        self.tableRelationships = self.getTableRelationships(tablesInfo)
        self.constraintTableConstraintKeyNames = getConstraintTablesConstraintKeyName(tablesInfo)

    def initialiseTables(self, tablesInfo):
        for t in tablesInfo:
            tableName = t['tableName']
            # the tables which contain links can have INTEGER columns, leave the rest as TEXT
            # (INTEGER might be faster?)
            if tableName == 'bgc_gcf' or tableName == 'mf_spectra'or tableName == 'spectra_bgc':
                columns = ', '.join(['{} INTEGER'.format(c) for c in t['tableData'][0].keys()])
            else:
                columns = ', '.join(['{} TEXT'.format(c) for c in t['tableData'][0].keys()])
                columns = []
                for c in t['tableData'][0].keys():
                    if c.endswith('_pk'):
                        columns.append('{} INTEGER PRIMARY KEY'.format(c))
                    else:
                        columns.append('{} TEXT'.format(c))
                columns = ', '.join(columns)

            sql = 'CREATE TABLE IF NOT EXISTS {} ({})'.format(t['tableName'], columns)
            self.db.execute(sql)
            if 'pk' in t['options']:
                sql = 'CREATE UNIQUE INDEX IF NOT EXISTS {}_index ON {} ({})'.format(t['tableName'], t['tableName'], t['options']['pk'])
                self.db.execute(sql)
            else:
                colnames = list(t['tableData'][0].keys())
                sql = 'CREATE UNIQUE INDEX IF NOT EXISTS {}_index on {} ({}, {})'.format(t['tableName'], t['tableName'], colnames[0], colnames[1])
                self.db.execute(sql)

        # only add data if needed
        cur = self.db.execute('SELECT COUNT({}) FROM {}'.format(tablesInfo[0]['options']['pk'], tablesInfo[0]['tableName']))
        rowcount = cur.fetchone()[0]
        if rowcount == 0:
            print('Adding data to tables')
            self.addNewData(tablesInfo)
        else:
            print('Tables already populated')

    def clearAlasqlTables(self, tablesInfo):
        for t in tablesInfo:
            self.db.execute('DELETE FROM {}'.format(t['tableName']))

    def addNewData(self, tablesInfo):
        for t in tablesInfo:
            # t['tablesData'] is a list of dicts, each representing a single row
            # {'col1': 'valcol1-1', 'col2': 'valcol2-1', ...},
            # {'col1': 'valcol1-2', 'col2': 'valcol2-2', ...}
            # etc
            cols = self.tableColumns[t['tableName']]
            with self.db:
                allvals = []
                stmt = 'INSERT INTO {} VALUES({})'.format(t['tableName'], ','.join('?' for c in cols))
                for row in t['tableData']:
                    allvals.append([row[c] for c in cols])
                self.db.executemany(stmt, allvals)
    
    def getFieldNames(self, tablesInfo):
        # this is meant to return a comma-separated list of '<table>.<col> AS table_col' 
        # entries for the 4 visible tables
        
        visible_tables = filter(isTableVisible, tablesInfo)
        firstdatarows = map(lambda t: {'tableName': t['tableName'], 'firstDataRow': t['tableData'][0]}, visible_tables)
        tabledata = map(lambda x: (x['tableName'], list(x['firstDataRow'].keys())), firstdatarows)
        all_cols = ['{}.{} AS {}_{}'.format(table, col, table, col) for table, cols in tabledata for col in cols]
        return ', '.join(all_cols)

    def assembleInnerJoinStatementFromRelationship(self, relationship, gen_to_met):
        def parseRelationship(r):
            if 'with' in r:
                if gen_to_met:
                    return 'INNER JOIN {} ON {}.{} = {}.{} '.format(r['tableName'], r['with'], r['using'], r['tableName'], r['using'])
                return 'INNER JOIN {} ON {}.{} = {}.{} '.format(r['with'], r['tableName'], r['using'], r['with'], r['using'])
            else:
                return ''

        rs = None
        if isinstance(relationship, list):
            rs = relationship
        else:
            rs = [relationship]

        innerJoinStatement = ''
        for i, r in enumerate(rs):
            innerJoinStatement += parseRelationship(r)

        return innerJoinStatement

    def makeSelectClause(self, tablesInfo, gen_to_met):
        fieldNames = self.getFieldNames(tablesInfo)
        if gen_to_met:
            return  'SELECT {} FROM {}'.format(fieldNames, tablesInfo[-1]['tableName'])
        return 'SELECT {} FROM {}'.format(fieldNames, self.firstTable)

    def makeInnerJoinClause(self, gen_to_met):
        if gen_to_met:
            return ' '.join(map(lambda x: self.assembleInnerJoinStatementFromRelationship(x, gen_to_met), reversed(self.tableRelationships)))
        return ' '.join(map(lambda x: self.assembleInnerJoinStatementFromRelationship(x, gen_to_met), self.tableRelationships))

    def getFirstTable(self, tablesInfo):
        return tablesInfo[0]['tableName']

    def getRelationship(self, tableInfo):
        def parseRelationship(r):
            return {'tableName': tableInfo['tableName'], 'with': r['with'], 'using': r['using']}

        if 'relationship' in tableInfo:
            if isinstance(tableInfo['relationship'], list):
                return list(map(parseRelationship, tableInfo['relationship']))
            else:
                return  {
                            'tableName': tableInfo['tableName'], 
                            'with': tableInfo['relationship']['with'], 
                            'using': tableInfo['relationship']['using']
                        }
        else:
            return {'tableName': tableInfo['tableName']}

    def getTableRelationships(self, tablesInfo):
        return list(map(self.getRelationship, tablesInfo))

    def getTableKeys(self):
        # construct a list of {'tableName': ..., 'tableKey': ...} dicts
        tks = []
        for t in self.tableRelationships:
            tk = {'tableName': t['tableName']}
            if 'using' in t:
                tk['tabkeKey'] = t['using']
            tks.append(tk)
        return tks

    def makeSQLquery(self, tablesInfo, skipConstraints, gen_to_met):
        selectClause = self.makeSelectClause(tablesInfo, gen_to_met)
        innerJoinClause = self.makeInnerJoinClause(gen_to_met)
        whereClause = self.makeWhereClause(tablesInfo, skipConstraints)
        return ' '.join([selectClause, innerJoinClause, whereClause])

    def makeWhereClause(self, tablesInfo, skipConstraints):
        whereSubClauses = self.makeWhereSubClauses()
        selectedWhereSubClauses = []
        for i, value in enumerate(whereSubClauses):
            if not skipConstraints[i]:
                selectedWhereSubClauses.append('{} IN @(?)'.format(whereSubClauses[i]))

        if len(selectedWhereSubClauses) > 0:
            whereSubClauses1 = ' AND '.join(selectedWhereSubClauses)
            return 'WHERE {}'.format(whereSubClauses1)
        else:
            return ''

    def makeWhereSubClauses(self):
        return list(map(lambda t: '{}.{}'.format(t['tableName'], t['constraintKeyName']), self.constraintTableConstraintKeyNames))

    def simpleQuery(self, tableName, tableFieldNames):
        sqlQuery = 'SELECT ' + ', '.join(tableFieldNames) + ' FROM {}'.format(tableName)
        resultset = self.db.execute(sqlQuery)

        data = []
        for row in resultset:
            data.append({f: row[f] for f in tableFieldNames})

        return data

    def prefixQuery(self, tableName, tableFieldNames):
        # need to add the prefix to each of the column names here...
        prefix = tableName + '_'
        tableFieldNames = list(map(lambda x: prefix + x, tableFieldNames))
        sqlQuery = 'SELECT DISTINCT ' + ', '.join(tableFieldNames) + ' FROM tempTable'

        resultset = self.db.execute(sqlQuery)

        # and then extract the results and remove the prefix from the names again
        data = []
        for row in resultset:
            data.append({f.replace(prefix, ''): row[f] for f in tableFieldNames})

        return data

    def createTempTable(self, tablesInfo, constraints):
        constraintTableNames = list(map(lambda t: t['tableName'], self.constraintTableConstraintKeyNames))
        unpackedConstraints = list(map(lambda n: constraints[n], constraintTableNames))
        skipConstraints = []
        selectedConstraints = []
        for i, uc in enumerate(unpackedConstraints):
            tableName = constraintTableNames[i]
            myTable = list(filter(lambda t: t['tableName'] == tableName, tablesInfo))[0]

            sc = False
            if len(uc) == 0 or len(uc) == len(myTable['tableData']):
                sc = True

            if not sc:
                selectedConstraints.append(uc)

            skipConstraints.append(sc)

        # HACK 1: check the query to see if we're going genomics->metabolomics or
        # metabolomics->genomics (so the JOINs can be performed in reverse order if needed)
        # skipConstraints is a 4-element list of bools corresponding to the 4 visible 
        # tables. if a value is true it means to ignore that table when constructing the
        # "constraints" for the query (the "WHERE" clauses). if BOTH of the first two
        # values are True, we're going met->gen, otherwise gen->met
        gen_to_met = skipConstraints[0] and skipConstraints[1]
        sqlQuery = self.makeSQLquery(tablesInfo, skipConstraints, gen_to_met)

        # HACK 2: replace any @(?) (alasql syntax) with standard placeholders
        c = 0
        while True:
            pos = sqlQuery.find('@(?)')
            if pos == -1:
                break
            # note final param, indicating to only replace 1 instance
            sqlQuery = sqlQuery.replace('@(?)', '({})'.format(','.join(['?' for i in range(len(selectedConstraints[c]))])), 1)
            c += 1

        # wrap this into a CREATE TEMP TABLE query
        finalQuery = 'CREATE TEMP TABLE tempTable AS {}'.format(sqlQuery)

        # flatten the nested list of constraints into a simple list so it can be substituted 
        # into the query as normal
        # TODO: note the casting to "int" here. for some reason the values can end
        # up being np.int64 instances and that breaks the lookup vs the integers in
        # the database (it gets treated as binary data). should figure out why that
        # is happening
        selectedConstraints = list(map(int, itertools.chain.from_iterable(selectedConstraints)))

        return self.db.execute(finalQuery, selectedConstraints)

    def deleteTempTable(self):
        return self.db.execute('DROP TABLE tempTable')

if __name__ == "__main__":
    # create some example data
    mf_df, mf_ds, mf_dt = get_data('mf', 5)
    spectra_df, spectra_ds, spectra_dt = get_data('spectra', 10)
    bgc_df, bgc_ds, bgc_dt = get_data('bgc', 10)
    gcf_df, gcf_ds, gcf_dt = get_data('gcf', 5)

    # create some example links
    # 0 is not used because it's for NA
    link_indices = [(1, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (4, 7), (5, 8), (5, 9), (5, 10)]
    mf_spectra = create_links(mf_df, spectra_df, 'mf_pk', 'spectra_pk', link_indices)

    link_indices = [(1, 1), (2, 2), (1, 3), (3, 4), (3, 5), (5, 7)]
    spectra_bgc = create_links(spectra_df, bgc_df, 'spectra_pk', 'bgc_pk', link_indices)

    link_indices = [(1, 1), (2, 1), (3, 1), (4, 1), (5, 2), (6, 3), (7, 4), (8, 5), (9, 5), (10, 5)]
    bgc_gcf = create_links(bgc_df, gcf_df, 'bgc_pk', 'gcf_pk', link_indices)

    # combine the data and link informations into a table info
    # this will be passed to the linker object when the Load button is clicked
    table_info = get_table_info(mf_df, mf_spectra, spectra_df, spectra_bgc, bgc_df, bgc_gcf, gcf_df)
    data_sources = {
        'mf_table': mf_ds,
        'spectra_table': spectra_ds,
        'bgc_table': bgc_ds,
        'gcf_table': gcf_ds
    }
    

    # sm = SqlManager(table_info)

    linker = Linker(table_info, 'test', 'test.sqlite3')


