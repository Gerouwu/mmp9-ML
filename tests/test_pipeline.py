"""Pruebas de umbrales, métricas, separación y preprocesamiento."""
import numpy as np
import pandas as pd
import pytest

from mmp9.data import activity_class, canonical_parent
from mmp9.experiment import metrics, preprocessor, split_indices, matrix, model, fit_model


def test_threshold_boundaries():
    assert [activity_class(x) for x in [1, 1000, 1001, 9999, 10000]] == [
        "active", "active", "intermediate", "intermediate", "inactive"]


def test_metrics_positive_is_active():
    result = metrics([0, 0, 0, 1, 1], [0, 0, 1, 1, 0])
    assert (result["tn"], result["fp"], result["fn"], result["tp"]) == (2, 1, 1, 1)
    assert result["accuracy"] == pytest.approx(0.6)
    assert result["f1"] == pytest.approx(0.5)
    assert result["sensitivity"] == pytest.approx(0.5)
    assert result["specificity"] == pytest.approx(2/3)


def test_no_negative_predictions_and_missing_class():
    result = metrics([0, 0, 1, 1], [1, 1, 1, 1])
    assert result["specificity"] == 0
    with pytest.raises(ValueError):
        metrics([1, 1], [1, 1])


def test_parent_smiles():
    assert canonical_parent("CCO.[Na+]") == canonical_parent("OCC") == "CCO"
    assert canonical_parent("invalid smiles") is None
    assert canonical_parent(None) is None


def test_all_splits_disjoint_reproducible():
    y = np.array([0] * 30 + [1] * 70)
    for seed in range(100):
        train, test, fit, valid = split_indices(y, seed)
        assert not set(train) & set(test)
        assert not set(fit) & set(valid)
        assert set(fit) | set(valid) == set(train)
        assert set(train) | set(test) == set(range(100))
        assert all(set(y[indices]) == {0, 1} for indices in [train, test, fit, valid])
        assert np.array_equal(test, split_indices(y, seed)[1])


def test_vocabulary_and_imputation_learn_only_train():
    vectorizer = preprocessor("SMILES")
    vectorizer.fit(pd.Series(["CCO", "CCO", "CCC", "CCC"]))
    before = vectorizer.vocabulary_.copy()
    vectorizer.transform(pd.Series(["ZZZZ"]))
    assert vectorizer.vocabulary_ == before and "ZZ" not in before
    prep = preprocessor("PaDEL")
    prep.fit([[1, 0], [3, 0], [np.nan, 0]])
    assert prep.named_steps["impute"].statistics_[0] == 2
    assert matrix(prep.transform([[1000, 0]])).shape == (1, 1)


@pytest.mark.parametrize("name", ["LogisticRegression", "RandomForest", "XGBoost"])
def test_models_produce_binary_probability(name):
    x = np.array([[0, 0], [0, 1], [1, 0], [1, 1]] * 10, dtype=np.float32)
    y = np.array([0, 0, 1, 1] * 10)
    estimator = fit_model(model(name, {}, 0), x, y, name)
    assert estimator.predict_proba(x).shape == (40, 2)
    assert set(estimator.classes_) == {0, 1}
