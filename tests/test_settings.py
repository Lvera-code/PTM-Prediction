from src.config.settings import Settings


def test_deepptmpred_threshold_for_tipo_calibrado():
    assert Settings.deepptmpred_threshold_for("phosphorylation") == 0.24020174


def test_deepptmpred_threshold_for_tipo_desconocido_usa_fallback():
    assert Settings.deepptmpred_threshold_for("tipo_inexistente") == Settings.DEEPPTMPRED_MIN_PROBABILITY


def test_deepptmpred_calibrated_thresholds_cubre_los_17_tipos():
    assert set(Settings.DEEPPTMPRED_CALIBRATED_THRESHOLDS) == set(Settings.DEEPPTMPRED_PTM_TYPES)
