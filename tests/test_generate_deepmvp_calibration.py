"""Tests de scripts/generate_deepmvp_calibration.py.

Solo cubre la logica pura (``trim_window``) y la orquestacion de descarga
(``download_all_data``, con red mockeada) -- ``generate_for_type`` importa
TensorFlow/el paquete ``lib`` de DeepMVP de forma diferida (dentro de la
funcion) y requiere el conda env dedicado ``deepmvp``, fuera de alcance para
la suite principal (mismo criterio que el resto de codigo dependiente de
motores externos, ver STATUS.md).
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from scripts.generate_deepmvp_calibration import MAX_WIDTH, download_all_data, trim_window


def test_trim_window_ancho_maximo_no_recorta():
    x = "A" * MAX_WIDTH
    assert trim_window(x, MAX_WIDTH) == x


def test_trim_window_recorta_centrado():
    # x siempre tiene exactamente MAX_WIDTH chars en la practica real (columna
    # 'x' de all_data.tar.gz, ver docstring del modulo) -- trim_window no
    # valida esto, asume el invariante.
    x = "".join(str(i % 10) for i in range(MAX_WIDTH))
    trimmed = trim_window(x, 11)
    trim = (MAX_WIDTH - 11) // 2
    assert trimmed == x[trim : MAX_WIDTH - trim]
    assert len(trimmed) == 11


@pytest.mark.parametrize("peptide_length", [9, 15, 21, 33, 51, 61])
def test_trim_window_longitud_de_salida_correcta(peptide_length):
    x = "M" * MAX_WIDTH
    assert len(trim_window(x, peptide_length)) == peptide_length


def test_download_all_data_no_descarga_si_ya_existe(tmp_path):
    dest_dir = tmp_path
    archive_path = dest_dir / "all_data.tar.gz"
    archive_path.write_bytes(b"fake")
    extracted = dest_dir / "all_data"
    extracted.mkdir()

    with patch("scripts.generate_deepmvp_calibration.urllib.request.urlretrieve") as mock_urlretrieve, \
         patch("scripts.generate_deepmvp_calibration.tarfile.open") as mock_tar_open:
        result = download_all_data(dest_dir)

    mock_urlretrieve.assert_not_called()
    mock_tar_open.assert_not_called()
    assert result == extracted


def test_download_all_data_descarga_y_extrae_si_falta(tmp_path):
    dest_dir = tmp_path / "work"

    def _fake_urlretrieve(url, filename):
        Path(filename).write_bytes(b"fake")

    with patch("scripts.generate_deepmvp_calibration.urllib.request.urlretrieve", side_effect=_fake_urlretrieve) as mock_urlretrieve, \
         patch("scripts.generate_deepmvp_calibration.tarfile.open") as mock_tar_open:
        result = download_all_data(dest_dir)

    mock_urlretrieve.assert_called_once()
    mock_tar_open.assert_called_once()
    assert result == dest_dir / "all_data"
    assert (dest_dir / "all_data.tar.gz").is_file()
