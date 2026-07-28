# PTM-Prediction

Pipeline de prediccion de zonas de modificacion post-traduccional (PTM) a
partir de FASTA o PDB/mmCIF. Segundo proyecto del CNB (encargado por Carlos
Oscar Sorzano), independiente de
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction).

## Arquitectura

Numeracion de fases alineada con el proyecto 1
([BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
Fase 1 -> Fase 1.5 (PDB unicamente) -> Fase 2 (motores) -> Fase 3 (nucleo).

- **Camino FASTA**: Fase 1 (saneamiento) -> Fase 2: **DeepMVP** (motor unico,
  6 tipos de PTM: fosforilacion, acetilacion, metilacion, sumoilacion,
  ubiquitinacion, N-glicosilacion).
- **Camino PDB**: Fase 1.5 (extraccion de secuencia ATMSEQ + mapeo de
  posiciones via `gemmi`) -> Fase 2: consenso **DeepMVP + DeepPTMPred** (17
  tipos de PTM; DeepPTMPred exige `pdb_path` obligatorio, sin modo
  solo-secuencia, de ahi la asimetria entre caminos).
- **Fase 3 (nucleo)**: anotacion/filtrado + logica de decision de flujo sobre
  las predicciones crudas de Fase 2 (`src/engines/ptm_annotation.py`).
  Fusiona consenso donde DeepMVP y DeepPTMPred coinciden en tipo+posicion.

## Estado actual (2026-07-27)

Pipeline completo end-to-end (Fase 1/1.5/3) implementado y con 67 tests.
Ningun motor tiene pesos instalados todavia en esta maquina — ver
`STATUS.md` para el detalle de que falta (pesos de DeepMVP, PyRosetta +
checkpoint ESM-2 de DeepPTMPred) y que ya se verifico instalable de verdad.

## Uso

```bash
python pipeline.py --input fasta_inputs/mi_proteina.fasta
python pipeline.py --input fasta_inputs/mi_estructura.pdb
```

## Instalacion

```bash
pip install -r requirements.txt
pytest tests/
```

DeepMVP y DeepPTMPred requieren venvs dedicados aparte (stacks
incompatibles entre si), nunca el venv principal del pipeline:

```bash
git clone https://github.com/bzhanglab/DeepMVP
conda create -n deepmvp python=3.7.10 -y
conda run -n deepmvp pip install -r DeepMVP/requirements.txt
export DEEPMVP_PYTHON_BIN=$(conda run -n deepmvp which python)
# Pesos: descarga manual desde https://deepmvp.ptmax.org/ -> DeepMVP/models/

git clone https://github.com/kuikui-wang/DeepPTMPred
# venv Python 3.10 (ver DeepPTMPred/pred/train_PTM/environment.yml),
# + pyrosetta-installer (ver README del repo) + checkpoint ESM-2 (~2.5GB)
```

Ver `STATUS.md` para el detalle completo de que ya se verifico (ambos
repos clonados y probados en esta maquina el 2026-07-27) y que falta.

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
