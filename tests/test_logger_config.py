"""Tests de src/utils/logger_config.py.

Cada test usa un nombre de logger unico (logging.getLogger es un singleton
global por nombre) para no interferir entre tests ni con los loggers reales
que ya crean otros modulos del pipeline al importarse.
"""

import logging

import src.utils.logger_config as logger_config


def test_setup_logger_es_idempotente_no_duplica_handlers(tmp_path, monkeypatch):
    monkeypatch.setattr(logger_config, "LOG_DIR", tmp_path / "logs")
    monkeypatch.setattr(logger_config, "LOG_FILE", tmp_path / "logs" / "ptm_pipeline.log")

    logger = logger_config.setup_logger("test_idempotente")
    n_handlers_primera_vez = len(logger.handlers)
    logger_again = logger_config.setup_logger("test_idempotente")

    assert logger is logger_again
    assert len(logger_again.handlers) == n_handlers_primera_vez


def test_console_handler_nivel_info(tmp_path, monkeypatch):
    monkeypatch.setattr(logger_config, "LOG_DIR", tmp_path / "logs")
    monkeypatch.setattr(logger_config, "LOG_FILE", tmp_path / "logs" / "ptm_pipeline.log")

    logger = logger_config.setup_logger("test_console_level")
    console_handlers = [h for h in logger.handlers if isinstance(h, logging.StreamHandler)
                        and not isinstance(h, logging.handlers.RotatingFileHandler)]

    assert len(console_handlers) == 1
    assert console_handlers[0].level == logging.INFO


def test_crea_directorio_de_log_y_handler_de_archivo(tmp_path, monkeypatch):
    log_dir = tmp_path / "logs_custom"
    monkeypatch.setattr(logger_config, "LOG_DIR", log_dir)
    monkeypatch.setattr(logger_config, "LOG_FILE", log_dir / "ptm_pipeline.log")

    logger = logger_config.setup_logger("test_crea_dir")

    assert log_dir.is_dir()
    file_handlers = [h for h in logger.handlers if isinstance(h, logging.handlers.RotatingFileHandler)]
    assert len(file_handlers) == 1


def test_logger_no_propaga_a_root(tmp_path, monkeypatch):
    monkeypatch.setattr(logger_config, "LOG_DIR", tmp_path / "logs")
    monkeypatch.setattr(logger_config, "LOG_FILE", tmp_path / "logs" / "ptm_pipeline.log")

    logger = logger_config.setup_logger("test_no_propaga")

    assert logger.propagate is False
    assert logger.level == logging.INFO
