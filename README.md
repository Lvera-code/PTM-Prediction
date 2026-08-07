# PTM-Prediction

[![tests](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml/badge.svg)](https://github.com/Lvera-code/PTM-Prediction/actions/workflows/tests.yml)

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
  `n_linked_glycosylation` esta deliberadamente EXCLUIDO de ESA fusion
  especifica (modelo de DeepPTMPred sin poder discriminativo real para este
  tipo, AUROC 0.51 -- ver STATUS.md), aunque ambos motores lo siguen
  reportando por separado -- en su lugar, Camino PDB tiene un consenso REAL
  distinto para este tipo especifico: **DeepMVP + EMNGly + StackGlyEmbed**
  (decision 2026-08-06, reemplaza al CoNglyPred original, sin pesos
  publicados -- ver STATUS.md "Decision 2"), `pasa_umbral` = al menos 1 de
  los motores disponibles pasa, `consenso` = al menos 2 de 3.
- **Corroboracion opcional** (nunca deciden `pasa_umbral`/`consenso`, solo
  informan sobre sitios que el consenso YA acepto): **MeToken** (tipo, Camino
  PDB), **StackGlyEmbed** (N-glicosilacion, Camino FASTA -- en Camino PDB es
  motor de consenso, ver arriba), **GlyGen** (evidencia experimental externa
  de glicosilacion).
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

## Estado actual (2026-08-06)

Pipeline completo end-to-end (Fase 1/1.5/2/3/A) implementado, instalado y
verificado con corridas REALES (no solo planeado) en esta maquina: DeepMVP
(pesos + calibracion real), DeepPTMPred (PyRosetta + ESM-2 + calibracion real
de los 17 tipos, 2 bugs de train/inferencia encontrados y corregidos --
phi/psi y plDDT, ver STATUS.md), MeToken, StackGlyEmbed y Fase A (ddG /
glicano / conjugacion, todos con corridas reales documentadas).

**Decision 2 (segundo motor de consenso para N-glicosilacion) CERRADA**:
CoNglyPred (candidato original) confirmado sin pesos publicados en ningun
sitio -- reemplazado por **EMNGly** (pesos reales, verificados a nivel de
bytes), implementado y wireado al consenso de `n_linked_glycosylation` junto
con StackGlyEmbed (promovido de corroboracion informativa a motor de
consenso en Camino PDB). Ver STATUS.md "Decision 2" para el detalle
completo, incluyendo los 2 go/no-go checks (ambos PASARON, 2026-08-07): el
alineamiento de `structure_emb` se verifico contra sitios reales de GlyGen
en un PDB con huecos (Alpha-1-antitrypsin, 1QLP), y el MCC publicado
(0.736) se reprodujo -- de hecho se supero, 0.8197 -- corriendo el pipeline
real sobre el set independiente de N-GlyDE. El umbral
`EMNGLY_MIN_PROBABILITY=0.5` ya no es provisional.

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

StackGlyEmbed (N-glicosilacion, AMBOS caminos -- FASTA y PDB). En Camino
FASTA sigue siendo puramente informativa (`Settings.STACKGLYEMBED_ENABLED`
degrada solo con un aviso si no esta disponible). En Camino PDB fue
PROMOVIDA 2026-08-06 a motor real de consenso junto con EMNGly (ver abajo),
reemplazando el rol que hubiera tenido DeepPTMPred si no estuviera
confirmado muerto para `n_linked_glycosylation` -- sigue siendo opcional
(el pipeline no se cae sin ella), pero deja de ser meramente informativa
ahi. A diferencia de los demas motores, NO requiere instalacion propia
aqui -- reusa el venv/pickles YA instalados en el proyecto hermano
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction)
(decision 2026-07-26: nunca se importa codigo entre proyectos, pero SI se
reusan recursos externos pesados ya instalados como venvs/pesos):

```bash
# Requiere que B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed
# ya exista (ver su propio README, Seccion 11). Si vive en otra ruta:
export STACKGLYEMBED_PYTHON_BIN=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed/bin/python
export STACKGLYEMBED_MODELS_DIR=/ruta/a/B-Cell-Epitope-Prediction/StackGlyEmbed/prediction
```

EMNGly (motor real de consenso para `n_linked_glycosylation`, Camino PDB
unicamente -- reemplaza a CoNglyPred, decision 2026-08-06, ver STATUS.md
"Decision 2": CoNglyPred confirmado sin pesos publicados en ningun sitio).
Opcional (`Settings.EMNGLY_ENABLED` degrada solo con un aviso, el consenso
de N-glicosilacion cae a DeepMVP+StackGlyEmbed si no esta disponible):

```bash
git clone https://github.com/StellaHxy/EMNgly
python3 -m venv .venv-emngly
.venv-emngly/bin/pip install fair-esm torch "scikit-learn==1.1.1" scipy pandas numpy tqdm wget
# 'wget' (paquete pip, no el binario CLI): dependencia transitiva real de
# MIF/sequence_models/trRosetta_utils.py (importado en cadena por
# pretrained.py aunque este proyecto nunca usa esa ruta de codigo) -- sin
# ella el import de MIF falla con ModuleNotFoundError (confirmado real
# 2026-08-07, no documentado en el environment.yml de EMNgly).
# pip resuelve numpy a la ultima 2.x por defecto, incompatible en runtime con
# el wheel compilado de scikit-learn==1.1.1 (ValueError: numpy.dtype size
# changed -- confirmado real 2026-08-07). Reinstalar despues, pinneado:
.venv-emngly/bin/pip install "numpy==1.23.5"
export EMNGLY_PYTHON_BIN=$(pwd)/.venv-emngly/bin/python

# MIF (embedding estructural, mif.pt) ya viene bundled en el clon -- sin
# descarga aparte.

# ESM-1b (~7.4GB) + companero de regresion de contactos (fair-esm SIEMPRE lo
# busca, mismo mecanismo que el checkpoint ESM-2 de DeepPTMPred arriba):
mkdir -p EMNgly/esm/checkpoints
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S.pt \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
curl -L -o EMNgly/esm/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt \
  https://dl.fbaipublicfiles.com/fair-esm/regression/esm1b_t33_650M_UR50S-contact-regression.pt

# SVM (verificado 2026-08-06 con HTTP Range real, ver STATUS.md): elegido
# N-GlyDE.pickle (negativos restringidos al sequon, el regimen correcto),
# no N-GlyAltas_classifier.pkl (misma carpeta, entrenado sobre el otro
# benchmark).
mkdir -p EMNgly/checkpoints
curl -L -o EMNgly/checkpoints/N-GlyDE.pickle \
  "https://drive.usercontent.google.com/download?id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t"
```

> Los pesos del SVM se generaron con `scikit-learn==1.1.1` -- el
> `environment.yml` del repo pinnea 1.5.1, pero cargar los pickles con una
> version distinta a la de entrenamiento dispara `InconsistentVersionWarning`
> (no garantizado). Pinnea `scikit-learn==1.1.1` en el venv, no lo que pide
> el `environment.yml` del repo.

Ver `STATUS.md` para el detalle completo de que ya se verifico (todos los
repos/envs clonados y probados en esta maquina, mas el recurso externo de
StackGlyEmbed, mas los 2 go/no-go checks de EMNGly -- ambos PASARON
2026-08-07) y que falta.

## Licencias

Este repositorio (codigo propio, ver `LICENSE`): CC BY-NC 4.0, uso no
comercial -- consistente con la dependencia real mas restrictiva (ver abajo).

- **DeepMVP**: GPL-3.0 (declarada en `DeepMVP/LICENSE` del repo original -- corregido
  2026-08-03, una version anterior de este README decia MIT por error, nunca
  verificado contra el archivo real hasta ahora).
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
- **EMNGly**: sin LICENSE declarado en ningun repo (`StellaHxy/EMNgly` ni su
  duplicado `Xiaoyang878/EMNgly`) -- el paper es CC-BY (Bioinformatics) pero
  eso NO cubre el codigo. Correo pendiente a Xiaoyang Hou/Shiwei Sun/Yaojun
  Wang (ICT-CAS), mismo patron que funciono con DeepPTMPred -- NO bloqueante
  (a diferencia de CoNglyPred, los pesos ya son publicos y descargables sin
  depender de esa respuesta). MIF (vendorizado dentro de EMNgly,
  `model/MIF/`) es de Microsoft (`microsoft/protein-sequence-models`),
  licencia BSD-2 permisiva verificada en el repo oficial -- la copia de
  EMNgly perdio su LICENSE al vendorizarlo.
- **PyRosetta** (Fase A): licencia academica/no-comercial de RosettaCommons,
  ya cubierta por el uso de investigacion/TFG de este proyecto.

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
