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

MeToken (corroboracion OPCIONAL e informativa del tipo, Camino PDB, ver
STATUS.md 2026-08-01 -- NUNCA un motor de consenso, se puede omitir sin
perder funcionalidad del pipeline principal, `Settings.METOKEN_ENABLED`
degrada solo con un aviso si no esta instalado):

```bash
git clone https://github.com/A4Bio/MeToken
conda create -n metoken python=3.10 -y
conda run -n metoken pip install --index-url https://download.pytorch.org/whl/cpu torch
conda run -n metoken pip install numpy scipy biopython transformers omegaconf tqdm pandas huggingface-hub h5py
conda run -n metoken pip install --no-build-isolation torch_scatter  # sin wheel prebuilt, compila desde fuente
export METOKEN_PYTHON_BIN=$(conda run -n metoken which python)
# Pesos: descarga release 1.0 (~88MB) y descomprime en MeToken/pretrained_model/
curl -L -o /tmp/pretrained_model.zip https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip
python -c "import zipfile; zipfile.ZipFile('/tmp/pretrained_model.zip').extractall('MeToken')"
```

Ver `STATUS.md` para el detalle completo de que ya se verifico (los 3 repos
clonados y probados en esta maquina) y que falta.

## Licencias

- **DeepMVP**: MIT (declarada en `DeepMVP/LICENSE` del repo original).
- **DeepPTMPred**: el repo (`github.com/kuikui-wang/DeepPTMPred`, vendorizado
  en `DeepPTMPred/`) no declara licencia propia, pero el paper asociado es
  CC BY-NC 4.0 (Oxford University Press). Se contacto directamente a los
  autores de correspondencia (Yong Liu, Junwen Wang) para confirmar que el
  codigo sigue los mismos terminos. **Junwen Wang confirmo por email el
  2026-07-29**: *"I confirm that the GitHub code follows the same CC BY-NC
  terms."* Uso no comercial (investigacion/TFG + integracion futura como
  plugin de Scipion por el CNB, institucion publica) encaja sin problema
  dentro de CC BY-NC. Detalle completo en `STATUS.md`.
- **MeToken**: MIT (declarada en `MeToken/LICENSE` del repo original).

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
