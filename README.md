# PTM-Prediction

Pipeline de prediccion de zonas de modificacion post-traduccional (PTM) a
partir de FASTA o PDB/mmCIF, con modelado estructural real (PyRosetta) de un
subconjunto de los sitios aceptados. Segundo proyecto del CNB (encargado por
Carlos Oscar Sorzano), independiente de
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction).

## Arquitectura

Numeracion de fases alineada con el proyecto 1
([BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)):
Fase 1 -> Fase 1.5 (PDB unicamente) -> Fase 2 (motores) -> Fase 3 (nucleo) -> Fase A
(modelado estructural, PDB unicamente).

- **Camino FASTA**: Fase 1 (saneamiento) -> Fase 2: **DeepMVP** (motor unico,
  6 tipos de PTM: fosforilacion, acetilacion, metilacion, sumoilacion,
  ubiquitinacion, N-glicosilacion) -> Fase 3.
- **Camino PDB**: Fase 1.5 (extraccion de secuencia ATMSEQ + mapeo de
  posiciones via `gemmi`) -> Fase 2: consenso **DeepMVP + DeepPTMPred** (17
  tipos de PTM; DeepPTMPred exige `pdb_path` obligatorio, sin modo
  solo-secuencia, de ahi la asimetria entre caminos) -> Fase 3 -> **Fase A**.
- **Fase 3 (nucleo)**: anotacion/filtrado + logica de decision de flujo sobre
  las predicciones crudas de Fase 2 (`src/engines/ptm_annotation.py`).
  Fusiona consenso donde DeepMVP y DeepPTMPred coinciden en tipo+posicion.
  `n_linked_glycosylation` esta deliberadamente EXCLUIDO del consenso (modelo
  de DeepPTMPred sin poder discriminativo real para este tipo, AUROC 0.51 --
  ver STATUS.md), aunque ambos motores lo siguen reportando por separado.
- **Corroboracion opcional** (nunca deciden `pasa_umbral`/`consenso`, solo
  informan sobre sitios que el consenso YA acepto): **MeToken** (tipo, Camino
  PDB), **StackGlyEmbed** (N-glicosilacion especificamente, ambos caminos,
  reusa el venv del proyecto hermano), **GlyGen** (evidencia experimental
  externa de glicosilacion).
- **Fase A (modelado estructural real via PyRosetta, Camino PDB unicamente)**:
  para un top-N de sitios por tipo (default 1, `Settings.FASE_A_TOP_N_PER_TYPE`
  -- modelar TODOS los sitios aceptados es computacionalmente inviable, un
  caso real como Tau acepta ~572), modela la estructura real segun 3 clases:
  parche quimico nativo + ddG (5 tipos: fosforilacion, acetilacion,
  hidroxilacion, gamma-carboxiglutamacion, metilacion de Lys), adjuncion +
  refinado de glicano (2 tipos: N/O-glicosilacion), o conjugacion isopeptidica
  real via `UBQ_GTPaseMover` (2 tipos: ubiquitinacion, sumoilacion). Los otros
  8/17 tipos no tienen modulo de Fase A (requeririan construir un residuo
  no-canonico propio, cheminformatica real, no implementado) -- se marcan
  `fase_a_estado="sin_soporte_fase_a"` en el reporte, nunca se omiten en
  silencio.

## Estado actual (2026-08-03)

Pipeline completo end-to-end (Fase 1/1.5/2/3/A) implementado, instalado y
verificado con corridas REALES (no solo planeado) en esta maquina: DeepMVP
(pesos + calibracion real), DeepPTMPred (PyRosetta + ESM-2 + calibracion real
de los 17 tipos, 2 bugs de train/inferencia encontrados y corregidos --
phi/psi y plDDT, ver STATUS.md), MeToken, StackGlyEmbed y Fase A (ddG /
glicano / conjugacion, todos con corridas reales documentadas). 154 tests
(`pytest tests/`).

Pendiente, no bloqueante: **CoNglyPred** como segundo motor de consenso
especifico para N-glicosilacion (candidato preferido, sin pesos publicados en
ningun sitio verificado -- correo enviado al autor de correspondencia,
StackGlyEmbed ya cubre parcialmente ese hueco mientras tanto). Ver STATUS.md,
seccion "Decision 2".

## Uso

```bash
python pipeline.py --input fasta_inputs/mi_proteina.fasta
python pipeline.py --input fasta_inputs/mi_estructura.pdb   # Camino PDB: motores + Fase A
```

## Instalacion

```bash
pip install -r requirements.txt
pytest tests/
```

DeepMVP, DeepPTMPred y Fase A requieren venvs/conda envs dedicados aparte
(stacks incompatibles entre si), nunca el venv principal del pipeline:

```bash
git clone https://github.com/bzhanglab/DeepMVP
conda create -n deepmvp python=3.7.10 -y
conda run -n deepmvp pip install -r DeepMVP/requirements.txt
export DEEPMVP_PYTHON_BIN=$(conda run -n deepmvp which python)
# Pesos: descarga manual desde https://deepmvp.ptmax.org/ -> DeepMVP/models/

git clone https://github.com/kuikui-wang/DeepPTMPred
conda create -n deepptmpred python=3.10 -y
# Instalar segun DeepPTMPred/pred/train_PTM/environment.yml (TensorFlow 2.15 +
# tensorflow-addons + torch + fair-esm) + pyrosetta-installer (ver README del
# repo) + checkpoint ESM-2 (~2.5GB, mas su companero *-contact-regression.pt)
export DEEPPTMPRED_PYTHON_BIN=$(conda run -n deepptmpred which python)

# Fase A (ddG / glicosilacion / ubiquitinacion-sumoilacion, src/structural/)
# REUSA el mismo conda env 'deepptmpred' -- ya tiene PyRosetta instalado, no
# hace falta un env adicional.
export FASE_A_PYTHON_BIN=$(conda run -n deepptmpred which python)
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

StackGlyEmbed (corroboracion OPCIONAL e informativa de N-glicosilacion,
AMBOS caminos -- FASTA y PDB, ver STATUS.md 2026-08-01 -- NUNCA un motor de
consenso, se puede omitir sin perder funcionalidad del pipeline principal,
`Settings.STACKGLYEMBED_ENABLED` degrada solo con un aviso si no esta
disponible). A diferencia de los demas motores, NO requiere instalacion
propia aqui -- reusa el venv/pickles YA instalados en el proyecto hermano
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)
(decision 2026-07-26: nunca se importa codigo entre proyectos, pero SI se
reusan recursos externos pesados ya instalados como venvs/pesos):

```bash
# Requiere que B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed
# ya exista (ver su propio README, Seccion 11). Si vive en otra ruta:
export STACKGLYEMBED_PYTHON_BIN=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed/bin/python
export STACKGLYEMBED_MODELS_DIR=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/prediction
```

Ver `STATUS.md` para el detalle completo de que ya se verifico (todos los
repos/envs clonados y probados en esta maquina, mas el recurso externo de
StackGlyEmbed) y que falta.

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
- **PyRosetta** (Fase A): licencia academica/no-comercial de RosettaCommons,
  ya cubierta por el uso de investigacion/TFG de este proyecto.

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
