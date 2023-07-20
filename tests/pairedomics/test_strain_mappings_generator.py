import json
from nplinker.pairedomics import extract_mappings_ms_filename_spectrum_id
from nplinker.pairedomics import \
    extract_mappings_original_genome_id_resolved_genome_id
from nplinker.pairedomics import extract_mappings_resolved_genome_id_bgc_id
from nplinker.pairedomics import extract_mappings_strain_id_ms_filename
from nplinker.pairedomics import extract_mappings_strain_id_original_genome_id
from nplinker.pairedomics import get_mappings_strain_id_bgc_id
from nplinker.pairedomics import get_mappings_strain_id_spectrum_id
from nplinker.pairedomics import podp_generate_strain_mappings
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain


def test_podp_generate_strain_mappings(monkeypatch, tmp_path):
    # mock functions called by the tested function
    mappings_strain_bgc = {
        "strain1": {"bgc1", "bgc2"},
        "strain2": {"bgc3"},
    }
    mappings_strain_spectrum = {
        "strain1": {"spec1", "spec2"},
        "strain2": {"spec3"}
    }

    monkeypatch.setattr(
        "nplinker.pairedomics.strain_mappings_generator.extract_mappings_strain_id_original_genome_id",
        lambda *args: {})  # any return value is fine
    monkeypatch.setattr(
        'nplinker.pairedomics.strain_mappings_generator.extract_mappings_original_genome_id_resolved_genome_id',
        lambda *args: {})
    monkeypatch.setattr(
        'nplinker.pairedomics.strain_mappings_generator.extract_mappings_resolved_genome_id_bgc_id',
        lambda *args: {})
    monkeypatch.setattr(
        "nplinker.pairedomics.strain_mappings_generator.get_mappings_strain_id_bgc_id",
        lambda *args: mappings_strain_bgc)

    monkeypatch.setattr(
        'nplinker.pairedomics.strain_mappings_generator.extract_mappings_strain_id_ms_filename',
        lambda *args: {})
    monkeypatch.setattr(
        'nplinker.pairedomics.strain_mappings_generator.extract_mappings_ms_filename_spectrum_id',
        lambda *args: {})
    monkeypatch.setattr(
        "nplinker.pairedomics.strain_mappings_generator.get_mappings_strain_id_spectrum_id",
        lambda *args: mappings_strain_spectrum)

    # Create the expected
    expected_dict = {
        "strain1": {"bgc1", "bgc2", "spec1", "spec2"},
        "strain2": {"bgc3", "spec3"}
    }
    expected_sc = StrainCollection()
    for strain_id, ids in expected_dict.items():
        strain = Strain(strain_id)
        for iid in ids:
            strain.add_alias(iid)
        expected_sc.add(strain)

    # Call function to generate strain mappings
    output_file = tmp_path / "output.json"
    result = podp_generate_strain_mappings("dummy_podp_project_file",
                                           "dummy_genome_status_file",
                                           "dummy_genome_bgc_mappings_file",
                                           "dummy_gnps_file_mapping_file",
                                           output_file)
    # check returned value
    assert isinstance(result, StrainCollection)
    assert result == expected_sc
    # check output file
    sc = StrainCollection.read_json(output_file)
    assert sc == expected_sc


def test_extract_mappings_strain_id_original_genome_id(tmp_path):
    test_data = {
        "genomes": [
            {
                "genome_label": "strain1",
                "genome_ID": {
                    "RefSeq_accession": "id1"
                }
            },
            {
                "genome_label": "strain1",
                "genome_ID": {
                    "RefSeq_accession": "id2"
                }
            },
            {
                "genome_label": "strain2",
                "genome_ID": {
                    "RefSeq_accession": "id3"
                }
            },
        ]
    }
    test_file = tmp_path / "test_data.json"
    with open(test_file, "w") as f:
        json.dump(test_data, f)

    expected_result = {
        "strain1": {"id1", "id2"},
        "strain2": {"id3"},
    }
    assert extract_mappings_strain_id_original_genome_id(
        test_file) == expected_result


def test_extract_mappings_original_genome_id_resolved_genome_id(tmp_path):
    test_data = {
        "genome_status": [
            {
                "original_id": "id1",
                "resolved_refseq_id": "refseq1",
                "resolve_attempted": True,
                "bgc_path": ""
            },
            {
                "original_id": "id2",
                "resolved_refseq_id": "refseq2",
                "resolve_attempted": True,
                "bgc_path": ""
            },
            {
                "original_id": "id3",
                "resolved_refseq_id": "refseq3",
                "resolve_attempted": True,
                "bgc_path": ""
            },
        ]
    }
    test_file = tmp_path / "test_data.json"
    with open(test_file, "w") as f:
        json.dump(test_data, f)

    expected_result = {"id1": "refseq1", "id2": "refseq2", "id3": "refseq3"}

    assert extract_mappings_original_genome_id_resolved_genome_id(
        test_file) == expected_result


def test_extract_mappings_resolved_genome_id_bgc_id(tmp_path):
    test_data = {
        "mappings": [
            {
                "genome_ID": "id1",
                "BGC_ID": ["bgc1", "bgc2"]
            },
            {
                "genome_ID": "id2",
                "BGC_ID": ["bgc3"]
            },
            {
                "genome_ID": "id3",
                "BGC_ID": []
            },
        ]
    }
    test_file = tmp_path / "test_data.json"
    with open(test_file, "w") as f:
        json.dump(test_data, f)
    expected_result = {
        "id1": {"bgc1", "bgc2"},
        "id2": {"bgc3"},
        "id3": set(),
    }
    assert extract_mappings_resolved_genome_id_bgc_id(
        test_file) == expected_result


def test_get_mappings_strain_id_bgc_id():
    # Test case 1: Test with empty mappings
    mappings_strain_id_original_genome_id = {}
    mappings_original_genome_id_resolved_genome_id = {}
    mappings_resolved_genome_id_bgc_id = {}
    expected_result = {}
    assert get_mappings_strain_id_bgc_id(
        mappings_strain_id_original_genome_id,
        mappings_original_genome_id_resolved_genome_id,
        mappings_resolved_genome_id_bgc_id) == expected_result

    # Test case 2: Test with one strain and one genome
    mappings_strain_id_original_genome_id = {"strain1": {"genome1"}}
    mappings_original_genome_id_resolved_genome_id = {
        "genome1": "resolved_genome1"
    }
    mappings_resolved_genome_id_bgc_id = {"resolved_genome1": {"bgc1"}}
    expected_result = {"strain1": {"bgc1"}}
    assert get_mappings_strain_id_bgc_id(
        mappings_strain_id_original_genome_id,
        mappings_original_genome_id_resolved_genome_id,
        mappings_resolved_genome_id_bgc_id) == expected_result

    # Test case 3: Test with multiple strains and genomes
    mappings_strain_id_original_genome_id = {
        "strain1": {"genome1", "genome2"},
        "strain2": {"genome3"},
        "strain3": {"genome4"}
    }
    mappings_original_genome_id_resolved_genome_id = {
        "genome1": "resolved_genome1",
        "genome2": "resolved_genome1",
        "genome3": "resolved_genome2",
        "genome4": ""
    }
    mappings_resolved_genome_id_bgc_id = {
        "resolved_genome1": {
            "bgc1",
        },
        "resolved_genome2": {"bgc2", "bgc3"},
    }
    expected_result = {"strain1": {"bgc1"}, "strain2": {"bgc2", "bgc3"}}
    assert get_mappings_strain_id_bgc_id(
        mappings_strain_id_original_genome_id,
        mappings_original_genome_id_resolved_genome_id,
        mappings_resolved_genome_id_bgc_id) == expected_result


def test_extract_mappings_strain_id_ms_filename(tmp_path):
    test_data = {
        "genome_metabolome_links": [
            {
                "genome_label": "strain1",
                "metabolomics_file": "http://example.com/file1.mzXML"
            },
            {
                "genome_label": "strain1",
                "metabolomics_file": "http://example.com/file2.mzXML"
            },
            {
                "genome_label": "strain2",
                "metabolomics_file": "http://example.com/file3.mzXML"
            },
            {
                "genome_label": "strain3",
                "metabolomics_file": "http://example.com/file4.mzXML"
            },
        ]
    }
    test_file = tmp_path / "test_data.json"
    with open(test_file, "w") as f:
        json.dump(test_data, f)
    expected_result = {
        "strain1": {"file1.mzXML", "file2.mzXML"},
        "strain2": {"file3.mzXML"},
        "strain3": {"file4.mzXML"},
    }

    assert extract_mappings_strain_id_ms_filename(test_file) == expected_result


def test_extract_mappings_ms_filename_spectrum_id(tmp_path):
    test_data = "cluster index\tAllFiles\nspec1\tfile1.mzXML:123###\nspec2\tfile2.mzXML:123###\nspec3\tfile2.mzXML:123###file3.mzXML:123###\n"
    test_file = tmp_path / "test_data.tsv"
    with open(test_file, "w") as f:
        f.write(test_data)
    expected_result = {
        "file1.mzXML": {"spec1"},
        "file2.mzXML": {"spec2", "spec3"},
        "file3.mzXML": {"spec3"}
    }

    assert extract_mappings_ms_filename_spectrum_id(
        test_file) == expected_result


def test_get_mappings_strain_id_spectrum_id():
    mappings_strain_id_ms_filename = {
        'strain1': {'file1.mzXML', 'file2.mzXML'},
        'strain2': {'file3.mzXML'},
        'strain3': {'file4.mzXML'}
    }
    mappings_ms_filename_spectrum_id = {
        'file1.mzXML': {'spec1'},
        'file2.mzXML': {'spec2', 'spec3'},
        'file3.mzXML': {'spec3'}
    }

    expected_mappings_dict = {
        'strain1': {'spec1', 'spec2', 'spec3'},
        'strain2': {'spec3'},
    }
    actual_mappings_dict = get_mappings_strain_id_spectrum_id(
        mappings_strain_id_ms_filename, mappings_ms_filename_spectrum_id)

    assert actual_mappings_dict == expected_mappings_dict
