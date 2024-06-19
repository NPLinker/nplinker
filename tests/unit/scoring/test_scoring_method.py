from nplinker.scoring import ScoringMethod


def test_has_value():
    assert ScoringMethod.has_value("metcalf") is True
    assert ScoringMethod.has_value("rosetta") is True
    assert ScoringMethod.has_value("nplclass") is True

    assert ScoringMethod.has_value("invalid") is False
