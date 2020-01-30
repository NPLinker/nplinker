class Linker {
    constructor(tablesInfo) {
        console.log('Linker initiliasing');
        console.log(tablesInfo);
        this.tablesInfo = tablesInfo;

        // precompute stuff that won't change
        this.defaultConstraints = this.getDefaultConstraints();
        this.fieldNames = this.getFieldNames();
        this.tableIdToIdColumnMap = this.getTableIdToIdColumnMap();

        console.log("IDCOL");
        console.log("%o", this.tableIdToIdColumnMap);

        // selections will change depending on what the user has chosen on screen
        this.selections = this.emptySelections();
        this.numSelected = this.countNumSelected();
        this.totalSelected = this.getTotalSelected()
        this.queryResult = null;
        this.sqlManager = new SqlManager(this.tablesInfo);
        console.log('Linker initiliased');
    }

    getDefaultConstraints() {
        return getConstraintTablesConstraintKeyName(this.tablesInfo)
            .reduce((constraints, tableInfo) => {
                constraints[tableInfo['tableName']] = this.getKeys(
                    this.tablesInfo, tableInfo['tableName'], tableInfo['constraintKeyName']);
                return constraints;
            }, {});
    }

    getFieldNames() {
        // Gets the field names for each visible table
        return this.tablesInfo
            .filter(isTableVisible)
            .map(tableInfo => ({
                'tableName': tableInfo['tableName'],
                'fieldNames': Object.keys(tableInfo['tableData'][0])
            }));
    }

    getTableIdToIdColumnMap() {
        // Get the table name and key used in the WHERE clause in the form tableName: key
        return getConstraintTablesConstraintKeyName(this.tablesInfo)
            .map(t => ({[t['tableName']]: t['constraintKeyName']}))
            .reduce((o, v) => Object.assign(o, v), {});
    }

    getKeys(tablesInfo, tableName, k) {
        // Gets the values of the key used in the table relationship for the SQL IN clause
        const data = tablesInfo
            .filter(isTableVisible)
            .filter(t => t['tableName'] === tableName)
            .map(t => t['tableData'])[0];

        const keys = data
            .map(d => d[k])
            .filter((k, idx, arr) => arr.indexOf(k) === idx);

        return keys;
    }

    getId(tableName, rowObject) {
        return rowObject[idColumn];
    }

    countNumSelected() {
        return getConstraintTablesConstraintKeyName(this.tablesInfo)
            .reduce((results, tableInfo) => {
                const tname = tableInfo['tableName'];
                results[tname] = this.selections[tname].length;
                return results;
            }, {});
    }

    getTotalSelected() {
        const values = Object.values(this.numSelected);
        const total = values.reduce((a, b) => a + b, 0);
        return total;
    }

    emptySelections() {
        return getConstraintTablesConstraintKeyName(this.tablesInfo)
            .reduce((results, tableInfo) => {
                const tname = tableInfo['tableName'];
                results[tname] = [];
                return results;
            }, {});
    }

    query() {
        // console.trace('query');

        // get alasql query result
        const constraints = this.makeConstraints();
        const resultset = this.sqlManager.queryDatabase(this.tablesInfo, constraints);

        // get results for each table
        const fieldNames = this.fieldNames;
        const tableData = {};
        for (let i = 0; i < fieldNames.length; i++) {
            const tableFieldNames = fieldNames[i];
            const tableName = tableFieldNames['tableName'];
            const data = this.sqlManager.prefixQuery(tableFieldNames, resultset);
            tableData[tableName] = data;
        }
        this.queryResult = tableData;
    }

    addConstraint(tableName, rowIndex, dataSource) {
        const idColumn = this.tableIdToIdColumnMap[tableName];
        const idVal = dataSource.data[idColumn][rowIndex]
        this.selections[tableName].push({
            idVal: idVal,
            rowIndex: rowIndex,
        });
        this.numSelected = this.countNumSelected();
        this.totalSelected = this.getTotalSelected();
    }

    removeConstraints(tableName) {
        this.selections[tableName] = [];
        this.numSelected = this.countNumSelected();
        this.totalSelected = this.getTotalSelected();
    }

    makeConstraints() {
        return getConstraintTablesConstraintKeyName(this.tablesInfo)
            .reduce((results, tableInfo) => {
                const tname = tableInfo['tableName'];
                results[tname] = this.selectionToConstraint(tname);
                return results;
            }, {});
    }

    selectionToConstraint(tableName) {
        if (this.numSelected[tableName] == 0) {
            return this.defaultConstraints[tableName];
        } else {
            return this.selections[tableName].map(x => x.idVal);
        }
    }

    updateDataSources(dataSources) {
        if (this.totalSelected > 0) {
            const fieldNames = this.fieldNames;
            const queryResult = this.queryResult;
            for (let i = 0; i < fieldNames.length; i++) { // update all the tables
                const tableFieldName = fieldNames[i];
                const tableName = tableFieldName['tableName'];
                if (this.numSelected[tableName] == 0) { // if table has no selection
                    // update table content
                    const dataSource = dataSources[tableName];
                    this.updateSingleTable(tableName, queryResult, dataSource);
                }
            }

        } else {
            this.resetTables(dataSources);
        }
    }

    updateSingleTable(tableName, queryResult, dataSource) {
        const newData = {};
        const tableQueryResult = queryResult[tableName]

        // natural sort in-place by id column
        // https://stackoverflow.com/questions/2802341/javascript-natural-sort-of-alphanumerical-strings
        const idColumn = this.tableIdToIdColumnMap[tableName];
        const collator = new Intl.Collator(undefined, {numeric: true, sensitivity: 'base'});
        tableQueryResult.sort((a, b) => collator.compare(a[idColumn], b[idColumn]));

        // convert from an array of objects to column data
        queryResult[tableName].forEach(function (row, i) {
            for (let prop in row) {
                if (row.hasOwnProperty(prop)) {
                    const value = row[prop];
                    if (!newData.hasOwnProperty(prop)) {
                        newData[prop] = []
                    }
                    newData[prop].push(value);
                }
            }
        })

        dataSource.data = newData;
        dataSource.change.emit();
    }

    resetTables(dataSources) {
        const fieldNames = this.fieldNames;
        const queryResult = this.queryResult;
        for (let i = 0; i < fieldNames.length; i++) { // update all the tables
            const tableFieldName = fieldNames[i];
            const tableName = tableFieldName['tableName'];
            const dataSource = dataSources[tableName];
            this.updateSingleTable(tableName, queryResult, dataSource);
        }
    }



}

const isTableVisible = tableInfo => tableInfo["options"]["visible"];

const getConstraintTablesConstraintKeyName = function(tablesInfo) {
    return tablesInfo
        .filter(isTableVisible)
        .map(t => ({'tableName': t['tableName'], 'constraintKeyName': t['options']['pk']}));
}

class SqlManager {

    constructor(tablesInfo) {
        this.initialiseAlasqlTables(tablesInfo);
        this.firstTable = this.getFirstTable(tablesInfo);
        this.tableRelationships = this.getTableRelationships(tablesInfo);
        this.constraintTableConstraintKeyNames = getConstraintTablesConstraintKeyName(tablesInfo);
    }

    initialiseAlasqlTables(tablesInfo) {
        tablesInfo.forEach(function (t) {
            // Create table
            let sql = "CREATE TABLE " + t['tableName'];
            console.log(sql);
            alasql(sql);
            // Create index
            if (t['options']['pk'] !== undefined) {
                sql = "CREATE UNIQUE INDEX tmp ON " + t['tableName'] + "(" + t['options']['pk'] + ")";
                console.log(sql);
                alasql(sql);
            }
            // Add data
            alasql.tables[t['tableName']].data = t['tableData'];
        });
    }

    clearAlasqlTables(tablesInfo) {
        tablesInfo.forEach(function (t) {
            alasql('DELETE FROM ' + t['tableName']);
        });
    }

    addNewData(tablesInfo) {
        tablesInfo.forEach(t => alasql.tables[t['tableName']].data = t['tableData']);
    }

    getFieldNames(tablesInfo) {
        return tablesInfo
            .filter(isTableVisible)
            .map(tableInfo => ({'tableName': tableInfo['tableName'], 'firstDataRow': tableInfo['tableData'][0]}))
            .map(tableData => Object.keys(tableData['firstDataRow'])
                .map(e => tableData['tableName'] + "." + e + " AS " + tableData['tableName'] + "_" + e))
            .reduce((fieldNamesArray, fieldNames) => fieldNamesArray.concat(fieldNames), [])
            .join(", ");
    }

    assembleInnerJoinStatementFromRelationship(relationship) {
        // debugger;
        function parseRelationship(r) {
            if (r['with']) {
                return "INNER JOIN " + r['with'] + " ON " + r['tableName'] + "." + r['using'] + " = " + r['with'] + "." + r['using'] + " ";
            } else {
                return "";
            }
        }

        let rs = undefined;
        if (relationship.constructor === Array) {
            rs = relationship; // an array of multiple relationships
        } else {
            rs = [relationship] // create an array of just one relationship
        }

        // process each relationship to make the final statement
        let innerJoinStatement = "";
        rs.forEach(function (r, i) {
            innerJoinStatement += parseRelationship(r);
        });
        return innerJoinStatement;
    }

    makeSelectClause(tablesInfo) {
        // Join each field into a select clause
        const fieldNames = this.getFieldNames(tablesInfo);
        // put the first table in the from clause
        const selectClause = "SELECT " + fieldNames + " FROM " + this.firstTable;
        return selectClause;
    }

    makeInnerJoinClause() {
        return this.tableRelationships
            .map(this.assembleInnerJoinStatementFromRelationship.bind(this))
            .join(" ");
    }

    getFirstTable(tablesInfo) {
        return tablesInfo[0]['tableName'];
    }

    getRelationship(tableInfo) {
        // debugger;
        function parseRelationship(r) {
            let parsed = {'tableName': tableInfo['tableName'], 'with': r['with'], 'using': r['using']};
            return parsed;
        }

        if (tableInfo['relationship']) {
            if (tableInfo['relationship'].constructor == Array) {
                // relationship is a list of dicts
                return tableInfo['relationship'].map(r => parseRelationship(r));
            } else {
                // relationship is a single dict
                return {
                    'tableName': tableInfo['tableName'],
                    'with': tableInfo['relationship']['with'],
                    'using': tableInfo['relationship']['using']
                }
            }
        } else {
            return {'tableName': tableInfo['tableName']};
        }
        ;
    }

    getTableRelationships(tablesInfo) {
        return tablesInfo
            .map(this.getRelationship);
    }

    getTableKeys() {
        // Returns the table name and the name of the key used in the where clause
        return this.tableRelationships
            .map(t => JSON.stringify({'tableName': t['tableName'], 'tableKey': t['using']}))
            .filter((tk, idx, tka) => tka.indexOf(tk) === idx)
            .map(t => JSON.parse(t));
    }

    makeSQLquery(tablesInfo, skipConstraints) {
        const selectClause = this.makeSelectClause(tablesInfo);
        const innerJoinClause = this.makeInnerJoinClause();
        const whereClause = this.makeWhereClause(tablesInfo, skipConstraints);
        return [selectClause, innerJoinClause, whereClause].join(" ");
    }

    makeWhereClause(tablesInfo, skipConstraints) {
        // add WHERE condition based on selected pks
        const whereSubClauses = this.makeWhereSubClauses();
        let selectedWhereSubClauses = [];
        whereSubClauses.forEach(function (value, i) {
            if (!skipConstraints[i]) {
                selectedWhereSubClauses.push(whereSubClauses[i] + ' IN @(?)');
            }
        });

        // combine whereSubClauses1 and whereSubClauses2
        if (selectedWhereSubClauses.length > 0) {
            const whereSubClauses1 = selectedWhereSubClauses.join(' AND ');
            return 'WHERE ' + whereSubClauses1;
        } else {
            return '';
        }
    }

    makeWhereSubClauses() {
        return this.constraintTableConstraintKeyNames
            .map(t => t['tableName'] + "." + t['constraintKeyName'])
    }

    queryDatabase(tablesInfo, constraints) {

        const constraintTableNames = this.constraintTableConstraintKeyNames.map(t => t['tableName']);
        const unpackedConstraints = constraintTableNames.map(n => constraints[n]);
        // console.log("unpackedConstraints.length = " + unpackedConstraints.length);
        let skipConstraints = [];
        let selectedConstraints = []
        unpackedConstraints.forEach(function (uc, i) {
            // find myTable matching by name
            let tableName = constraintTableNames[i];
            let myTable = tablesInfo.filter(t => t['tableName'] === tableName)[0];

            // if the where subClause includes ALL the data of that table, then skip it
            let sc = false
            if (uc.length == 0 || uc.length == myTable['tableData'].length) {
                sc = true;
            }
            if (!sc) {
                selectedConstraints.push(uc);
            }
            skipConstraints.push(sc);
            // console.log('%d. skip %s (%s)', i, sc, uc);
        });

        const sqlQuery = this.makeSQLquery(tablesInfo, skipConstraints);
        console.log(sqlQuery);
        const compiledSQLQuery = alasql.compile(sqlQuery);

        return compiledSQLQuery(selectedConstraints);
    }

    prefixQuery(tableFieldNames, dataSource) {
        const tableName = tableFieldNames['tableName'];
        const prefix = tableName + '_';
        const fieldNames = tableFieldNames['fieldNames'].map(x => prefix + x);

        const sqlQuery = "SELECT DISTINCT " + fieldNames.join(", ") + " FROM ?";
        console.log(sqlQuery);
        const temp = alasql(sqlQuery, [dataSource]);

        temp.map(x => { // for each row in the sql results
            Object.keys(x).map(key => { // rename the properties to remove the table name in front
                const newkey = key.replace(prefix, '');
                x[newkey] = x[key];
                delete (x[key]);
            });
        });

        return temp;
    }

}
