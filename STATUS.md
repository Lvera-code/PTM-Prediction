# STATUS

Estado actual del proyecto, no un diario de sesiones — reescribir en vez de
acumular. Ver `01-Proyectos/PTM-Prediction/` en el vault para el historial
de decisiones.

## DECISIÓN ARQUITECTÓNICA — CERRADA 2026-07-27

**Confirmado por Enzo: se mantiene la arquitectura actual, sin cambios de
código.** DeepMVP sigue siendo el único motor de Camino FASTA; DeepPTMPred
sigue siendo el motor de Camino PDB (estructura real); Camino PDB sigue
corriendo DeepMVP+DeepPTMPred en consenso, sin sumar un tercer motor por
ahora. MTPrompt-PTM se evaluó a fondo y se descarta como reemplazo — motivo
completo abajo. Nada en `pipeline.py`/`src/engines/` necesitaba cambiar:
esto ya era el diseño implementado, la decisión era si mantenerlo o no.

Enzo pidió reconsiderar que motor cubre Camino FASTA y cual Camino PDB,
evaluando **MTPrompt-PTM** (`github.com/hanye311/MTPrompt-PTM`) como
candidato. Verificado con el mismo nivel de profundidad que DeepMVP/
DeepPTMPred el 2026-07-27 (repo clonado real, venv `mtprompt` creado,
`pip install -r requirements.txt` sin errores, CLI real confirmado:
`python test.py --config_path ... --model_path ... --data_path <fasta>
--PTM_type <tipo> --save_path ...`).

**Hallazgos clave:**
- MTPrompt-PTM es **solo-secuencia** (input `protein_sequence`, confirmado
  en `config/PTM_config_prompt_tuning_test.yaml`) — NO usa el PDB/
  coordenadas 3D en absoluto durante la inferencia, pese a llamarse
  "Structure-Aware": su "conciencia de estructura" viene de como se
  PREENTRENÓ el backbone (S-PLM v2), no de features estructurales que se
  calculen por proteina en tiempo de inferencia (a diferencia de
  DeepPTMPred, que SI calcula SASA/phi/psi/secundaria reales via
  PyRosetta sobre el PDB de cada proteina). Consecuencia: si sustituye a
  DeepPTMPred, el pipeline pierde el ÚNICO motor que aprovecha estructura
  3D real -- Camino PDB pasaria a ser, en la práctica, dos motores de
  secuencia en consenso, no secuencia+estructura.
- No requiere PyRosetta (confirmado, ausente de `environment.yml`/
  `requirements.txt`) — resuelve el bloqueante real de hoy.
- Instala limpio (Python 3.11, ~3 min, sin conflictos reales — un solo
  warning cosmético de version de `packaging`).
- **Hallazgo real nuevo, no trivial**: el propio código de MTPrompt-PTM
  (`model.py::prepare_adapter_h_model`) construye su encoder ESM-2 via
  `esm_adapterH.pretrained.esm2_t33_650M_UR50D()`, que internamente llama
  `load_model_and_alphabet_HUB` (no la rama local que si usa nuestro
  runner de DeepPTMPred) — **por defecto SI intenta descargar de red**
  (`dl.fbaipublicfiles.com`) la primera vez que corre, a menos que el
  checkpoint ya este cacheado en `~/.cache/torch/hub/checkpoints/` de
  antemano. Parcheable (mismo patron que ya se aplico para DeepPTMPred:
  forzar la rama local pre-poblando el cache o pasando ruta `.pt`
  directamente), pero NO es "100% local de fabrica" como se asumio al
  proponerlo -- requeriria el mismo tipo de ajuste que ya se hizo para
  DeepPTMPred.
- Pesos propios (`best_model_13ptm_final.pth`, un unico archivo
  multi-tarea): descarga manual desde Google Drive, no descargados aun.
  README tambien ofrece Docker (`hanye0311/mtprompt:v1`) como alternativa
  de despliegue, no probado.
- Cubre 13 tipos (le faltan los 6 que DeepPTMPred si tiene: hydroxylation,
  gamma_carboxyglutamic_acid, glutarylation, glutathionylation,
  s_nitrosylation, citrullination; suma palmitoilacion, que ni DeepMVP ni
  DeepPTMPred cubren). Licencia MIT (mejor que DeepPTMPred, sin licencia
  declarada).

**Decisión final (2026-07-27): opción (a), sin cambios.** Motivo: DeepMVP
tiene evidencia verificada mas solida (benchmark real contra 8 herramientas
en Nature Methods) y cubre exactamente los 6 tipos ya decididos el 26/07 —
los 7 tipos extra de MTPrompt-PTM quedan fuera de ese alcance, no son una
ventaja real para este proyecto. DeepPTMPred sigue siendo el unico motor
con estructura 3D real, algo que se perderia si se sustituyera. No se
rediseño `ptm_annotation.py` porque no hubo cambio de motores.

MTPrompt-PTM y MusiteDeep (alternativas FASTA) y LkaM-PTM/SAPP
(alternativas PDB) quedan documentadas como candidatos evaluados y
descartados, no lineas de trabajo activas — ver el desglose completo en
`01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-implementacion-fase1-a-fase3-nucleo.md`
del vault y el artifact publicado en la sesion del 27-07 (comparacion de
arquitecturas de los 6 motores investigados).

## Hecho — pipeline completo end-to-end (2026-07-27), motores sin instalar

Las 3 fases estan implementadas y conectadas en `pipeline.py`. Ningun motor
(DeepMVP, DeepPTMPred) esta instalado todavia en esta maquina — correr el
pipeline hoy falla con un error accionable en Fase 3 (`DeepMVPExecutionError`
apenas invoca DeepMVP, confirmado con una corrida real). 67 tests
(`pytest tests/`), sin binarios/modelos externos, subprocess mockeado donde
aplica.

- **Fase 1** (`src/utils/fasta_parser.py`): saneamiento FASTA. Politica de
  residuos no canonicos RELAJADA (decision de Enzo 2026-07-27) para igualar
  la tolerancia real de DeepMVP (verificada en `lib/PeptideEncode.py`): ya
  no rechaza fatal ningun caracter, solo avisa para los que DeepMVP degrada
  (fuera de 20 estandar + U + B; 'X' no cuenta, DeepMVP lo trata como "sin
  señal"). A proposito distinto de `BCell-Epitope-Prediction` (rechazo
  fatal ahi, correcto porque BepiPred-3.0 SI aborta): cada Fase 1 refleja
  la tolerancia real de su motor.
- **Fase 1.5** (`src/utils/structure_parser.py`): ATMSEQ + mapeo de
  posiciones desde PDB/mmCIF via `gemmi`. Identico a proyecto 1 (logica
  generica).
- **Enrutador de input** (`src/utils/input_router.py`): identico a proyecto 1.
- **DeepMVPEngine** (`src/engines/deepmvp_engine.py`, motor unico Camino
  FASTA / motor 1 de 2 Camino PDB): subprocess sobre
  `github.com/bzhanglab/DeepMVP`, verificado leyendo el repo directamente
  (README.md, `DeepMVP.py`, `lib/PTModels.py`, `lib/Metrics.py`). CLI real
  con flags propios (`predict -m <model_dir> -d <fasta> -t 2 -o <out_dir>`),
  salida fija `site_prediction.tsv`.
- **DeepPTMPredEngine** (`src/engines/deepptmpred_engine.py` +
  `src/engines/_deepptmpred_runner.py`, motor 2/2 del consenso, Camino PDB
  unicamente): verificado leyendo `github.com/kuikui-wang/DeepPTMPred`
  directamente (README.md, `predict.py`, `e2_single_data.py`,
  `environment.yml`). A diferencia de DeepMVP, ESTE REPO NO TIENE CLI —
  `ptm_type`/`pdb_path`/`protein_id`/ruta del checkpoint ESM estan
  hardcodeados dentro de su bloque `if __name__ == "__main__":`. Se
  construyo un runner propio (`_deepptmpred_runner.py`) que importa las
  clases parametrizadas de `predict.py` (`PredictConfig`, `PTMPredictor`) y
  REIMPLEMENTA la extraccion de features ESM-2 en vez de usar
  `e2_single_data.py::extract_full_sequence_esm`, que tiene un bug real
  confirmado: redefine `custom_checkpoint_path` como variable LOCAL con una
  ruta absoluta hardcodeada de AutoDL (`/root/autodl-tmp/...`), ignorando
  cualquier valor pasado. Predice UN tipo de PTM por invocacion (17 en
  total), no todos a la vez como DeepMVP — el engine invoca 17 veces por
  accession y concatena.
- **Nucleo de Fase 3** (`src/engines/ptm_annotation.py`, B: anotacion +
  D: filtro): implementado segun el diseno de
  `01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md`.
  Correspondencia de tipos DeepMVP↔DeepPTMPred verificada leyendo ambos
  repos (no asumida): de los 6 tipos biologicos de DeepMVP (8 modelos por
  residuo), 7 tienen equivalente en DeepPTMPred; **`phosphorylation_y`
  (fosforilacion en Y) NO tiene ningun modelo equivalente en DeepPTMPred**
  (su unico modelo de fosforilacion cubre S/T, no Y) — hallazgo nuevo, ni
  siquiera en Camino PDB hay consenso posible para ese tipo. Los 10 tipos
  exclusivos de DeepPTMPred (no ~11 como se estimaba el 26: hydroxylation,
  gamma_carboxyglutamic_acid, malonylation, crotonylation, succinylation,
  glutathionylation, s_nitrosylation, glutarylation, citrullination,
  o_linked_glycosylation) se incluyen en el nucleo, marcados
  `consenso=false`. Umbral: `fpr <= DEEPMVP_MAX_FPR` para DeepMVP,
  `probability >= DEEPPTMPRED_MIN_PROBABILITY` (0.5, PROVISIONAL — DeepPTMPred
  no expone ningun mecanismo de calibracion, a diferencia del `fpr` de
  DeepMVP) para DeepPTMPred-solo; en filas con consenso, `pasa_umbral` es la
  UNION (pasa si al menos uno de los dos pasa su propio umbral), `consenso`
  es la INTERSECCION (True solo si ambos pasan).
- **Orquestador** (`pipeline.py`): las 3 fases conectadas end-to-end en
  ambos caminos. Camino FASTA maneja FASTA multi-accession correctamente
  (anota cada accession con su propia secuencia, no las concatena).
  Confirmado con corrida real: Fase 1 se completa y persiste
  (`<accession>_clean.fasta`), Fase 3 falla con `DeepMVPExecutionError`
  accionable (repo no clonado) — no falla en silencio a mitad de camino.
- Repo local (`git init`), sin remoto todavia — se crea en GitHub
  (`Lvera-code/PTM-Prediction`, publico) cuando haya algo funcional
  end-to-end CON los motores reales instalados.

## Instalacion real — verificado 2026-07-27 (repos clonados, dependencias instaladas de verdad)

Ambos repos se clonaron de verdad en esta maquina (`git clone`, dentro de
`PTM-Prediction/`, gitignorados) para verificar friccion de instalacion
real, no solo leer codigo. Contenido identico al verificado por API/curl
horas antes (`diff` byte a byte contra `predict.py` confirmado).

- **DeepMVP — INSTALADO, CON PESOS, VERIFICADO CON UNA PREDICCION REAL
  (2026-07-28).** Entorno conda real `deepmvp`
  (`/home/enzo/miniconda3/envs/deepmvp/bin/python`, Python 3.7.10):
  `pip install -r requirements.txt` completo sin errores (TensorFlow
  2.4.2 CPU, numpy 1.19.5, todo el stack). `python DeepMVP.py predict -h`
  imprime exactamente el CLI documentado arriba. Verificado tambien el
  caso de error real: con `-m` apuntando a una carpeta de modelos VACIA,
  DeepMVP.py revienta con `ValueError: No objects to concatenate` (traza
  cruda de pandas, opaca) — confirma que el chequeo proactivo de
  `DeepMVPEngine._validate_installation` (falla antes, con mensaje
  accionable) es necesario, no cosmetico. Pesos reales (`models.tar.gz`
  del Shiny app, NO `all_data.tar.gz` — ese es solo el dataset de
  train/test) instalados en `DeepMVP/models/`, los 8 tipos completos
  (10 modelos .h5 por tipo, ensemble). **Corrida real end-to-end sobre el
  FASTA de ejemplo del propio repo** (`example/Q5S007.fasta`):
  `site_prediction.tsv` final con 1292 filas reales, scores sensatos para
  los 8 tipos. Un mensaje "The file for computing FPR doesn't exist" sale
  una vez por tipo en la primera corrida (recalcula el umbral porque no
  habia un `site_prediction.tsv` previo dentro de `models/<tipo>/`) — no
  bloqueante, confirmado corriendo dos veces. Warning benigno de numpy/API
  mismatch tambien presente, tampoco afecta la ejecucion.
- **DeepPTMPred — entorno conda real construido y verificado, ESM-2
  funcional de punta a punta (2026-07-28).** Entorno conda `deepptmpred`
  (`/home/enzo/miniconda3/envs/deepptmpred`, Python 3.10), instalado desde
  `pred/train_PTM/environment.yml` (sin `cudatoolkit`/`cudnn`, esta maquina
  no tiene GPU): `tensorflow==2.15.0` + `tensorflow-addons==0.23.0` +
  `torch==2.13.0` + `fair-esm==2.0.0` instalan e importan sin error (mismo
  `UserWarning` de deprecacion de TFA ya conocido, no bloqueante).
  Checkpoint ESM-2 (`esm2_t33_650M_UR50D.pt`, 2.6GB, descarga manual del
  usuario) y su companero `esm2_t33_650M_UR50D-contact-regression.pt`
  (3.7KB) instalados en `DeepPTMPred/esm/checkpoints/`; ambos verificados
  con `torch.load()` (claves `args`/`cfg`/`model` presentes, no es una
  descarga corrupta). **`_extract_esm_features` de
  `_deepptmpred_runner.py` ejecutado de verdad** (no mockeado) contra el
  checkpoint real con una secuencia de prueba de 167 residuos: devuelve un
  array `(167, 1280)` — longitud de secuencia × dimension de embedding de
  `esm2_t33_650M`, exactamente lo esperado. **Unico bloqueante real
  restante para DeepPTMPred: PyRosetta.** Se intento la instalacion
  automatica real via `pip install pyrosetta-installer` +
  `install_pyrosetta()` (el metodo que documenta el propio README de
  DeepPTMPred) — el paquete instalador SI se instala, pero la descarga del
  wheel real falla en esta maquina: el mirror por defecto
  (`west.rosettacommons.org`) responde `404` en la ruta que el instalador
  consulta (`.../latest.html`), y el mirror alternativo (`graylab.jhu.edu`)
  falla la verificacion TLS (cadena de certificados no confiable en este
  entorno). Pendiente: usuario descargando el wheel manualmente
  (`PyRosetta4.Release.python310.ubuntu.wheel`, build mas reciente de
  2026-07-25) para instalar con `pip install <wheel>` sin depender de esos
  mirrors.
- **RESUELTO 2026-07-29: licencia de DeepPTMPred confirmada por el propio
  autor de correspondencia.** El repo no tenia LICENSE declarado, pero el
  paper si es CC BY-NC 4.0 (Oxford University Press) -- se le pregunto
  directamente a los autores de correspondencia (Yong Liu, Junwen Wang) si
  el codigo del repo sigue los mismos terminos. **Junwen Wang respondio
  2026-07-29 confirmando que si**: "I confirm that the GitHub code follows
  the same CC BY-NC terms." Uso no comercial (investigacion/TFG, e
  integracion futura como plugin de Scipion por el CNB, institucion
  publica) encaja dentro de CC BY-NC sin problema. Sin mas pasos
  pendientes en este punto.
- **Confirmado 100% local en toda la cadena (2026-07-27, re-verificado con
  ejecucion real 2026-07-28)**: verificado leyendo `esm/pretrained.py` de
  github.com/facebookresearch/esm directamente y ahora tambien con una
  corrida real. `pretrained.load_model_and_alphabet(path)` solo llama a
  red (`dl.fbaipublicfiles.com`) si el argumento NO termina en `.pt` (rama
  hub); como el runner siempre pasa una ruta `.pt` local, siempre entra por
  la rama `load_model_and_alphabet_local` (`torch.load()` puro sobre
  disco) — confirmado en la corrida real de arriba, sin acceso a red.
- El runner propio (`_deepptmpred_runner.py`) tiene su pieza de extraccion
  de features ESM-2 verificada contra el entorno real (ver arriba). Sigue
  sin poder correrse de punta a punta (prediccion completa via
  `PredictConfig`/`PTMPredictor`) porque falta PyRosetta -- solo esa pieza
  sigue probada unicamente con `subprocess.run` mockeado en tests.

## PyRosetta instalado y DeepPTMPred 100% funcional end-to-end (2026-07-28, tarde)

Usuario descargo el wheel manualmente (`pyrosetta-2026.30+release.bc091c65b8-cp310-cp310-linux_x86_64.whl`,
1.5GB) a `/mnt/c/Users/USUARIO/Downloads/` (WSL) y confirmo la instalacion.
Instalado con `pip install <wheel>` en el conda env `deepptmpred`,
verificado con `pyrosetta.init()` real (banner + version impresos, sin
error). Bloqueante que quedaba abierto desde el 27-07, cerrado.

**Corrida real end-to-end de `_deepptmpred_runner.py`** (no mockeada,
`AF-P10636-F1-model_v4.pdb` = Tau/MAPT, tipo `phosphorylation`): revelo y
arreglo 3 bugs reales adicionales, ninguno anticipado antes de tener
PyRosetta real corriendo:
- **Dependencias de `environment.yml` incompletas en el conda env real**:
  `matplotlib`, `seaborn`, `scikit-learn`, `imbalanced-learn`, `tqdm`,
  `joblib`, `logomaker` faltaban pese a estar listadas en
  `pred/train_PTM/environment.yml` (el `pip install -r requirements.txt`
  original del 27-07 no las trajo). Instaladas manualmente, confirmado
  `import predict` ya no falla.
- **Bug real en `predict.py::PTMPredictor.__init__`** (nueva clase de bug,
  distinta a la de `e2_single_data.py` ya documentada): el modelo guardado
  tiene una capa `Lambda` (`model.py::182`,
  `Lambda(lambda xin: K.sum(xin, axis=1))`) cuya funcion serializada
  referencia el simbolo `K` (alias de `tensorflow.keras.backend`) — Keras
  SOLO resuelve esos simbolos via el diccionario `custom_objects` pasado a
  `load_model`, nunca via los globals del modulo que la importa (aunque
  `predict.py` si tiene `K` en su propio namespace). El `custom_objects`
  real de `PTMPredictor.__init__` no incluye `'K'` -> `NameError: name 'K'
  is not defined` al cargar cualquier modelo. Confirmado real corriendo sin
  parche (revienta) y con parche (carga bien) antes de aplicarlo. Parcheado
  en `src/engines/_deepptmpred_runner.py::_load_predict_module` (envuelve
  `predict.load_model`, inyecta `'K': predict.K` en `custom_objects`) — NO
  se edita `predict.py` (vendored), mismo criterio que el resto del runner.
- Con ambos arreglos: corrida real completa, PyRosetta calcula SASA real +
  extrae plDDT del B-factor, ESM-2 carga desde cache local, el modelo
  predice — 130 filas de salida real para
  `phosphorylation` sobre Tau, scores 0.13-0.91, sensatos.

**`pipeline.py` completo (Camino PDB, consenso) corrido de verdad** sobre
el mismo caso: descubrio un bloqueante REAL adicional, distinto y anterior
al de DeepPTMPred, en el lado de DeepMVP -- ver seccion siguiente.

## RESUELTO: DeepMVP no calculaba `fpr` — calibracion generada de verdad (2026-07-28, noche)

Bloqueante encontrado horas antes en esta misma sesion (ver commit
`9fddb8c`): `DeepMVPEngine`/el nucleo de Fase 3 dependen de la columna
`fpr`, que `DeepMVP.py` solo calcula si existe
`DeepMVP/models/<tipo>/site_prediction.tsv` (archivo de VALIDACION con
columna `y` real) -- ausente en `models.tar.gz` para los 8 tipos,
confirmado ademas como bug conocido upstream (issues #1/#2 de
`bzhanglab/DeepMVP` en GitHub, dos usuarios externos independientes
reportaron lo mismo).

**Enzo pidio verificar mas a fondo antes de lanzar una descarga masiva de
UniProt** ("¿buscaste exhaustivamente? los datos deberían estar en DeepMVP,
¿no?"). Al revisar con mas cuidado, la respuesta era si: **no hacia falta
ningun FASTA externo**. La columna `x` de `all_data.tar.gz` (encontrado
antes, ver commit `9fddb8c`) ya viene con el ancho MAXIMO que pide
cualquier submodelo del ensemble (61 = flank 30, verificado leyendo los 8
`model.json`: ningun submodelo de ningun tipo pide `peptide_length` > 61).
Los submodelos que piden una ventana mas angosta se derivan recortando esa
misma ventana centrada (mismo criterio que
`DeepMVP/lib/DataIO.py::getPeptideSequence`) -- el requisito de un FASTA
`db` completo solo aparece si se usa el CLI (`DeepMVP.py predict -i`, que
llama incondicionalmente a `processing_prediction_data`, otro bug real
confirmado ahi), no si se reusa directamente la logica interna de
codificacion/ensemble (`lib.PeptideEncode.encodePeptides`,
`lib.Utils.combine_rts`).

Construido `scripts/generate_deepmvp_calibration.py` (standalone, conda
env `deepmvp`): descarga `all_data.tar.gz`, para cada uno de los 8 tipos
corre los `.h5` del ensemble ya instalado sobre `<tipo>_testing_70.tsv`
(recortando la ventana por submodelo como se describe arriba), promedia
con `combine_rts` (identico al `ptm_predict` real) y escribe
`DeepMVP/models/<tipo>/site_prediction.tsv`. **Verificado con AUROC real
por tipo, no solo "corre sin error"**: acetylation_k 0.986,
glycosylation_n 0.996, methylation_k 0.906, methylation_r 0.932,
phosphorylation_st 0.991, phosphorylation_y 0.992, sumoylation_k 0.974,
ubiquitination_k 0.989 -- coincide con las cifras publicadas en el paper,
confirma que la calibracion generada es correcta.

**Pipeline completo (Camino PDB) corrido de verdad end-to-end con la
calibracion ya en su sitio**: `AF-P10636-F1-model_v4.pdb` (Tau) -> 572
sitios PTM reales pasan el umbral (17 tipos representados), 98 con
consenso real DeepMVP+DeepPTMPred. Primera vez que el pipeline completa de
punta a punta sin fallar. `Camino FASTA` (mismo filtro `fpr`) queda
tambien desbloqueado por el mismo fix.

Sin resolver (menor, no bloqueante): que subconjunto usar de los 3
disponibles (`testing_70/80/90`, sufijos no documentados en el codigo del
repo -- muy probablemente umbrales de identidad de secuencia estilo CD-HIT
usados en el paper para medir generalizacion a distinta redundancia con el
training set, INFERIDO por convencion estandar en este tipo de papers, no
confirmado leyendo el texto del paper). El script usa `testing_70` (el mas
conservador, menos redundante con training) por defecto, parametrizable
via `--testing-suffix`.

## Fase A clase 1 + utilidad compartida + Extension 3 (ΔΔG) — implementadas y verificadas 2026-07-28

Construidas en el orden decidido el 28-07 (utilidad compartida -> Fase A
clase 1 -> Extension 3), en `src/structural/` (standalone, requiere
`pyrosetta`, nunca se importa desde el paquete `src` principal — mismo
patron que `_deepptmpred_runner.py`).

**Hallazgo real que corrige la estimacion "12/17 tipos" de la decision del
28-07** (esa cifra asumia cobertura por analogia con fosforilacion, sin
inspeccionar el enum/patches reales todavia — la propia nota lo dejaba como
pendiente explicito): verificado por DOS vias independientes contra el
PyRosetta ya instalado —
1. `ls database/chemical/residue_type_sets/fa_standard/patches/*.txt` (real).
2. Introspeccion en caliente de `pyrosetta.rosetta.core.chemical.VariantType`.

Solo **5 de los 17 tipos** tienen un `VariantType`/patch nativo listo para
usar: `phosphorylation` (S/T/Y), `acetylation` (K), `hydroxylation` (P),
`gamma_carboxyglutamic_acid` (E), `lys_methylation` (K, solo mono/di/tri
sin distinguir grado -- los motores no lo distinguen). Los otros 8 tipos de
"clase 1" (malonylation, arg_methylation, crotonylation, succinylation,
glutathionylation, s_nitrosylation, glutarylation, citrullination) NO
tienen patch nativo (grep de esos terminos sobre el directorio completo de
patches: cero resultados) -- requeririan construir un residuo no-canonico
propio (params file con topologia/cargas parciales/icoor), tarea de
cheminformatica real, NO implementada, NO fabricada.

- `src/structural/pyrosetta_ptm_patch.py`: `apply_ptm_patch` (aplica el
  variant real via `add_variant_type_to_pose_residue`, valida que el
  residuo real coincida con el esperado por tipo) + `relax_neighborhood`
  (FastRelax cartesiano, 1 repeat, restringido a un radio del sitio via
  `NeighborhoodResidueSelector` + `MoveMapFactory` + `PreventRepackingRLT`
  sobre el resto). CLI que aplica un parche real y opcionalmente relaja,
  volcando un PDB de salida.
- Verificado con una ejecucion real: acetilacion sobre LYS24 de
  `AF-P10636-F1-model_v4.pdb` -> `LYS` pasa a `LYS:acetylated`
  (`pose.residue(24).name()` antes/despues).
- `src/structural/ddg_estimate.py` (Extension 3): reutiliza la utilidad de
  arriba, compara `ref2015_cart` de la pose WT relajada localmente contra
  la pose parcheada relajada localmente (mismo protocolo en ambas).

  **Robustez de produccion — RESUELTO 2026-07-28 noche**: `FastRelax` es
  estocastico (Monte Carlo), asi que una sola relajacion puede caer en un
  minimo local peor de lo real. Se anadio `nstruct` (default 3, practica
  estandar de `cartesian_ddg`/`ddg_monomer` de Rosetta): corre N
  trayectorias INDEPENDIENTES por estado (cada una desde una pose fresca,
  no encadenadas) y se queda con la de menor energia por estado -- nunca
  promedia scores directamente, eso mezclaria trayectorias mal convergidas
  con las buenas. Tambien reporta la desviacion estandar entre trayectorias
  como medida de confianza.

  Verificado con una corrida real (`nstruct=2`, mismo sitio de acetilacion
  en LYS24 de Tau): ambas trayectorias convergieron al MISMO score exacto
  (1901.6247898992478 WT, 1845.4015068481324 mutado, std=0.0) -- confirmado
  que el generador de numeros aleatorios de PyRosetta si esta bien
  sembrado por proceso (semillas distintas en corridas separadas,
  verificado con `pyrosetta.rosetta.numeric.random.rg().get_seed()`), asi
  que esta convergencia identica no es un bug de semilla fija: es una
  propiedad real de este caso (vecindario pequeno, pocos residuos movibles,
  paisaje energetico simple) mas que evidencia de robustez general -- otros
  sitios/radios mayores podrian mostrar mas varianza entre trayectorias.
  Mecanismo verificado correcto (minimo entre trayectorias reales
  independientes, no una ilusion de "N corridas" que en realidad sean la
  misma). **Salvedad que sigue en pie**: Tau es intrinsecamente desordenada
  (plDDT medio 49.34 en este modelo AlphaFold) -- ningun numero de
  repeticiones arregla que la region de partida tenga baja confianza
  estructural, solo reduce el ruido de muestreo del relax en si.
- Ninguno de los dos scripts esta conectado a `pipeline.py` todavia --
  sigue la decision 2026-07-27 de que D (el filtro) no rutea a Extension
  3/Fase A porque esas fases no eran parte del nucleo. Se invocan por
  separado, manualmente o desde un futuro orquestador de Fase A/Extension 3.

## Fase A clase 2 (glicosilacion) — IMPLEMENTADA Y VERIFICADA (2026-07-28, noche)

`GlycanTreeModeler` REFINA la conformacion de un glicano ya adjunto -- para
adjuntarlo desde cero hacia falta otro Mover, `SimpleGlycosylateMover`
(mismo namespace `protocols.carbohydrates`, tambien confirmado real por
introspeccion). Ninguno de los dos motores (DeepMVP/DeepPTMPred) predice
la COMPOSICION del glicano, solo el sitio -- se usa el nucleo biosintetico
conservado como default documentado (no una prediccion):
`N-glycan_core` (Man3GlcNAc2, presente en TODO N-glicano maduro, unico
default universalmente defendible) para N-linked, `core_1_O-glycan`
(antigeno T, el mas comun mamiferos pero no universal -- existen 8 cores
distintos en `common_glycans/`) para O-linked. Ambos strings IUPAC
verificados leyendo los `.iupac` reales del PyRosetta instalado.

`src/structural/pyrosetta_glycan_patch.py`: `attach_glycan` (valida
residuo real N/S/T antes de parchear) + `refine_glycan`
(`GlycanTreeModeler`, rounds configurable). Verificado con dos corridas
reales sobre `AF-P10636-F1-model_v4.pdb`: N-glycan_core en ASN484 (sequon
real confirmado NAT) -- 5 residuos de azucar anadidos, `ASN484` pasa a
`ASN:N-glycosylated`; core_1_O-glycan en THR52 -- 2 residuos anadidos,
`THR52` pasa a `THR:O-conjugated`. Refinamiento con `GlycanTreeModeler`
(1 round) corre sin error, ~64s para el caso N-glycan_core.

## Fase A clase 3 (ubiquitinacion/sumoilacion) — IMPLEMENTADA (2026-08-01), corrige la conclusion de 2026-07-28

Enzo pidio terminar todo lo posible. Se investigo el mecanismo real de
Rosetta para conjugacion covalente de proteina completa (isopeptidico
Lys-NZ + Gly-Cterm de ubiquitina/SUMO), no solo se descarto por analogia:

- **Lado de la Lisina: SI existe y funciona.** Patch `SidechainConjugation`
  (real, cargado por defecto, confirmado con introspeccion de
  `ResidueTypeSet.patches()`) convierte LYS en
  `LYS:sidechain_conjugation` -- elimina los 3 hidrogenos de NZ y abre un
  punto de conexion real (`ADD_CONNECT NZ`), listo para enlazar a otra
  cadena.
- **Lado de la ubiquitina/SUMO: NO disponible en este build, con motivo
  documentado por los propios desarrolladores de Rosetta.** Existe un
  patch construido especificamente para esto --
  `patches/branching/C-terminal_conjugation.txt`, cuyo propio comentario
  dice literalmente *"currently only implemented for glycine conjugated
  ubiquitin"* -- pero esta **comentado/desactivado** en el manifiesto real
  `chemical/residue_type_sets/fa_standard/patches.txt` (linea 85):
  `#patches/branching/C-terminal_conjugation.txt  # Something is broken
  with this patch. ~Labonte` (Labonte = colaborador real de RosettaCommons
  especializado en quimica de carbohidratos). Confirmado con introspeccion
  real (`ResidueTypeSet.patches()` sobre una pose con ubiquitina real
  cargada, `1UBQ.pdb` de RCSB): el patch NO esta entre los 122 patches
  cargados por defecto, mientras que `SidechainConjugation` si lo esta.
- Esto CORROBORA de forma independiente y mas fuerte el hallazgo de la
  investigacion del 28-07 por la manana (foro de RosettaCommons: "poor fit
  for the isopeptide domain") -- no es solo que Rosetta tenga baja
  precision ahi, es que la unica herramienta nativa construida
  especificamente para esto esta desactivada por sus propios
  desarrolladores por estar rota, sin indicar que especificamente falla.
- **Decision: no se intento reconstruir el mecanismo a mano** (cirugia de
  pose de bajo nivel: eliminar OXT manualmente, declarar el enlace via
  `Conformation.declare_chemical_bond`, idealizar geometria sin ninguna
  referencia para validar que el resultado es quimicamente correcto). Sin
  el patch oficial como referencia de la geometria/bookkeeping de valencia
  correcta, y sin capacidad de validar el resultado contra una estructura
  real de conjugado (no existe superposicion trivial disponible), el
  riesgo de producir una estructura con apariencia plausible pero
  quimicamente incorrecta (bond lengths/angles no idealizados
  correctamente, o un problema de valencia no detectado) es alto y no
  detectable sin revision experta -- inaceptable para un pipeline de
  investigacion real. Ubiquitina real (`1UBQ.pdb`, descargada de RCSB,
  76 residuos, Gly75-Gly76 en el C-terminal) queda disponible en
  `DeepPTMPred/data/` (o donde se guarde) para cuando se decida seguir por
  esta via con mas tiempo, posiblemente cargando el patch manualmente
  (`-chemical:patch_selectors` apuntando al archivo directamente pese a
  estar excluido del manifiesto) y validando el resultado con cuidado.

67 tests (`pytest tests/`) seguian pasando en ese momento (2026-07-28) --
ningun test nuevo para `src/structural/` (requiere `pyrosetta`, ausente del
entorno donde corre la suite principal, mismo criterio que el resto de
codigo dependiente de motores externos: verificacion real documentada aqui
en vez de mock). **Actualizado 2026-07-29: 97 tests** tras cerrar los items
1/2/8 de la auditoria de robustez mas el cliente GlyGen (ver seccion de
auditoria abajo y "Proximos pasos" item 3) -- `src/structural/*.py` que
dependen de `pyrosetta` siguen sin tests, mismo motivo; `glygen_client.py`
es la unica excepcion (no requiere `pyrosetta`, si tiene tests).

### CORRECCION 2026-08-01: la conclusion de arriba era incorrecta -- IMPLEMENTADA via `UBQ_GTPaseMover`

La investigacion del 28-07 solo miro el patch (`C-terminal_conjugation.txt`,
desactivado por Rosetta) y concluyo "no hay via nativa disponible en este
build". Un subagente Opus (investigacion de robustez/produccion, pedida por
Enzo) encontro leyendo el .cc real de
`RosettaCommons/rosetta` (`source/src/protocols/chemically_conjugated_docking/UBQ_GTPaseMover.cc`
y la app `UBQ_Gp_LYX-Cterm.cc`) que existe una via COMPLETAMENTE
INDEPENDIENTE del patch roto: `UBQ_GTPaseMover`, un `protocols::moves::Mover`
normal que PyRosetta expone automaticamente (cero compilacion), que usa un
`ResidueType` completo (`LYX`) en vez del patch para el enlace isopeptidico.
Verificado real, no solo leido:

- `pyrosetta.rosetta.protocols.chemically_conjugated_docking.UBQ_GTPaseMover`
  se instancia sin error contra el wheel YA instalado (`deepptmpred`, PyRosetta
  2026.30) -- cero instalacion nueva.
- Requiere DOS flags de init no obvias, ambas confirmadas por prueba y error
  real contra el mismo error que reporta Rosetta si falta cualquiera de las
  dos: `-chemical:exclude_patches SidechainConjugation` (la misma flag que usa
  el test de integracion oficial de Rosetta para esta app) + `-extra_res_fa
  <ruta a LYX.params>` (sin esta segunda, `ResidueTypeSet` revienta con "The
  residue LYX could not be generated" pese a tener la primera flag puesta).
- Corrida real end-to-end contra `DeepPTMPred/data/AF-P10636-F1-model_v4.pdb`
  (Tau, LYS24, mismo residuo que ya usa Extension 3 para el ddG de
  acetilacion): `initialize()` solo -> pose de 834 residuos (758 Tau + 76
  ubiquitina), residuo objetivo mutado a `LYX`. `apply()` completo
  (`refine_cycles=3`, valor rapido de desarrollo) -> `MS_SUCCESS` en 16.4s,
  5 metricas de torsion de la cola movil recalculadas con valores reales no
  triviales. Mismo mecanismo verificado tambien para SUMO1 (pose de 857 =
  758 + 99 residuos) -- un solo Mover sirve para ambos tipos, solo cambia el
  PDB de referencia.
- Dos hallazgos reales de integracion (documentados en el docstring del
  modulo, no solo aqui): (1) `initialize()` escribe siempre
  `starting_pose_EMPTY_JOB_use_jd2_0000.pdb` en el cwd actual del proceso
  (confirmado, no un supuesto) porque `JobDistributor::current_job()` fuera
  de una corrida real via JobDistributor devuelve un job dummy en vez de
  `nullptr` -- no revienta, pero hay que confinarlo a un directorio temporal
  desechable; (2) las 5 metricas por-decoy que Rosetta escribe en ese mismo
  job dummy son inaccesibles desde Python -- se recalculan directamente
  desde la pose de salida via `mover.get_atom_IDs()` (expuesto, real) en vez
  de intentar leerlas del job.
- Un tercer hallazgo real, no de integracion sino de seguridad de tipo:
  Rosetta mismo NO valida que el residuo objetivo sea una Lisina antes de
  forzarlo a `LYX` (el `runtime_assert` real para eso esta comentado en el
  .cc) -- validado en nuestro wrapper en su lugar. Tambien se valida que la
  pose sea de una sola cadena (`runtime_assert` real de Rosetta, la
  violacion revienta alto en vez de silenciosa) y que la numeracion PDB
  coincida 1:1 con la numeracion de pose (`initialize()` hace una doble
  conversion pdb2pose sobre un indice que ya es de pose -- solo correcta si
  ambas numeraciones coinciden; verificado que es el caso para los PDBs de
  AlphaFold de una sola cadena que usa este pipeline, pero no asumido para
  cualquier PDB futuro).

Implementado en `src/structural/ubiquitin_sumo.py` (ver su docstring para el
detalle completo verificado linea por linea del .cc real). PDBs de
referencia empaquetados en `src/structural/data/`: `ubiquitin_reference.pdb`
(RCSB 1UBQ, el mismo que usa el propio equipo de Rosetta en sus tests de
integracion) y `sumo1_reference.pdb` (RCSB 1A5R, modelo 1 de la NMR,
recortado programaticamente -- buscando el motivo real "QTGG" en la
secuencia, no un numero de residuo asumido -- hasta la diglicina C-terminal
madura real, descartando 4 residuos de cola de maduracion sin procesar que
tiene el constructo original). Sirve tambien para SUMO sin cambios de
codigo (mismo mecanismo Gly-Gly terminal). Sin tests nuevos en la suite
principal (mismo motivo que el resto de `src/structural/`: requiere
`pyrosetta`; 101 tests siguen pasando sin cambios, verificado). Vault:
decision pendiente de escribir en `01-Proyectos/PTM-Prediction/Decisiones/`.

## Calibracion real de DeepPTMPred + bug de distribucion phi/psi (2026-07-30/31)

**2026-07-30**: implementado `scripts/generate_deepptmpred_calibration.py`
(corre el modelo real -- los `.h5` ya instalados via `predict.PTMPredictor`,
el mismo codigo que `_deepptmpred_runner.py` usa en produccion -- sobre una
muestra real de `DeepPTMPred/data/ptm_data.csv.gz`, secuencia+estructura real
por proteina descargada de AlphaFold DB). AUROC/umbral (Youden's J) medidos
por cada uno de los 17 tipos soportados, no asumidos. `scripts/
_rebuild_calibration_summary.py` agrega los 17 TSV por tipo en un
`summary.tsv` sin re-correr el modelo (necesario porque el script principal
solo escribe `summary.tsv` una vez, al final de su propia corrida -- dos
corridas parciales se pisaban entre si).

**2026-07-31: bug real encontrado y corregido, afecta la fiabilidad de TODO
lo calibrado hasta ese momento.** Los AUROC medidos salian sistematicamente
muy por debajo de los publicados en el paper (Briefings in Bioinformatics,
DOI bbag321) para practicamente los 17 tipos -- 2 de ellos (`hydroxylation`
0.34, `lys_methylation` 0.46) incluso por debajo de azar. Verificado leyendo
`DeepPTMPred/pred/train_PTM/data_loader.py` linea 139-141: `phi_center`/
`psi_center` se calculan con `half_window = (max(window_sizes)-1)//2 = 25`,
pero el array `phi`/`psi` del CSV de origen solo tiene 11 elementos -- la
condicion `len(x) > half_window` es SIEMPRE falsa, asi que el modelo
shippeado nunca vio nada distinto de 0.0 en esos dos angulos, ni en
entrenamiento ni en el test set que produce los AUC del paper (`trainer.py`
usa el mismo `load_dataset()` para train Y test). Pero `predict.py`'s
`PyRosettaCalculator.calculate_features` (inferencia REAL, la que usa
`_deepptmpred_runner.py` en produccion) SI calcula phi/psi reales via
PyRosetta -- un mismatch real train/inferencia, no un bug de esta calibracion.

Verificado empiricamente forzando `phi=psi=0.0` en inferencia (monkeypatch
de `PyRosettaCalculator.calculate_features`, aplicado en
`_deepptmpred_runner.py::_load_predict_module`, el mismo punto ya usado
para el parche real de `'K'`/`custom_objects`): AUROC de `hydroxylation`
sube 0.342->0.934 (paper: 0.965), `lys_methylation` 0.462->0.883 (paper:
0.899) -- recupera el rendimiento publicado. SASA/estructura
secundaria/plDDT NO estan afectados (solo phi_center/psi_center), asi que
el modelo no queda ciego a estructura, solo a la geometria fina de phi/psi
-- un techo de precision real y permanente de estos pesos concretos
(irreversible sin reentrenar con el bug de `data_loader.py` corregido, fuera
de alcance de este proyecto).

**Todos los datos de calibracion anteriores a este fix son invalidos**
(median el comportamiento roto). Movidos a
`DeepPTMPred/data/calibration_STALE_prephifix_2026-07-31/` (gitignored, no
borrados). Relanzada la calibracion completa de los 17 tipos a n=75 con el
parche puesto -- **terminada 2026-07-31**, resultado en
`DeepPTMPred/data/calibration/summary.tsv` (gitignored):

```
ptm_type                     auroc     suggested_threshold  n_proteins
gamma_carboxyglutamic_acid   0.970133  0.214867              44
hydroxylation                0.934372  0.375260              69
malonylation                 0.888000  0.422340             145
lys_methylation              0.882629  0.410402             143
glutathionylation            0.872889  0.472519             145
ubiquitination                0.870222  0.543494             149
phosphorylation              0.859970  0.243480             137
arg_methylation               0.850450  0.343200             139
sumoylation                   0.840180  0.409133             139
acetylation                   0.833538  0.629997             140
o_linked_glycosylation        0.823447  0.310599             121
succinylation                 0.808711  0.482645             144
s_nitrosylation                0.683243  0.500678             139
glutarylation                  0.673067  0.461801             103
citrullination                 0.657401  0.391692              46
crotonylation                  0.606957  0.935762              34
n_linked_glycosylation         0.490667  0.997414             110
```

14/17 tipos quedan en rango solido-a-excelente, muchos igualando o
superando el paper (confirma que el bug de phi/psi era la causa real de
los AUROC bajos, no motores debiles). 4 tipos quedan mediocres
(s_nitrosylation/glutarylation/citrullination/crotonylation, 0.61-0.68 --
`crotonylation` con muestra pequena, solo 34 proteinas unicas). **Un tipo
sigue realmente roto: `n_linked_glycosylation` (AUROC 0.49, peor que
azar) incluso tras el fix** -- no investigado todavia, causa distinta al
bug de phi/psi (que ya se corrigio), pendiente de mirar aparte, no
bloqueante para el resto.

### 2026-08-01: investigacion de `n_linked_glycosylation` (AUROC 0.49) y los 4 tipos mediocres (via subagente Opus)

**`n_linked_glycosylation`: CONFIRMADO limite real de los pesos publicados,
no un bug nuestro.** Tres confirmaciones independientes:

1. `DeepPTMPred/data/calibration/n_linked_glycosylation_calibration.tsv`
   esta saturado, no invertido: los 125 sitios puntuan 0.913-0.999,
   `prediction=1` en TODAS las filas -- el modelo no tiene salida
   discriminativa alguna (descarta un bug de label flip).
2. El propio repo trae sus metricas de entrenamiento para este tipo
   (`DeepPTMPred/pred/train_PTM/result/results_n_linked_glycosylation_esm2_kfold/ptm_data_200_29_64_average_results.txt`):
   **AUC 0.4950**, FP Rate 0.94, matriz de confusion de un fold sin NINGUN
   verdadero negativo (`[[0,5],[0,51]]`). Nuestra medicion (0.4907)
   reproduce esto casi exacto.
3. Ni siquiera se recupera aplicando el fix de plDDT de abajo (0.4907 ->
   0.5064, sigue siendo azar).

Causa raiz probable, verificada contra `ptm_data.csv.gz`: este tipo tiene
2115 positivos vs 355 negativos (85.6% positivo, el mas desbalanceado de
los 17 junto a crotonylation, 84.1%) -- y **el SMOTE del propio
`trainer.py:236` es un no-op para los 17 tipos** (`target_samples =
min(2000, class_counts[1])`, el `min` siempre da el conteo actual, nunca
sobremuestrea de verdad; contradice lo que dice el Methods del paper). El
paper mismo reporta 0.616 para este tipo (Table 1) y lo atribuye a
"limited negative sample diversity within conserved motifs" -- los
negativos son asparraginas DENTRO de sequones N-X-[S/T] intactos, el mismo
motivo de secuencia que verian los positivos. Los pesos shippeados ni
siquiera alcanzan ese 0.616 -- reproducen el artefacto 0.495.

**Recomendacion del subagente**: excluir `n_linked_glycosylation` del
consenso de DeepPTMPred (no es un umbral mal calibrado, es un modelo sin
poder discriminativo). El umbral de 0.997 en `summary.tsv` para este tipo
no significa nada. Refuerza la opcion ya listada como item 5 abajo
(cross-validacion con StackGlyEmbed del proyecto 1 para este tipo
especifico).

**4 tipos mediocres: NO es un techo del modelo, es un bug real de
inferencia adicional, mismo patron que el de phi/psi.** `predict.py:260`
calcula `local_plDDT` iterando POR ATOMO
(`atom.get_bfactor() for atom in structure.get_atoms()`), pero
`predict.py:299-303` indexa esa lista POR NUMERO DE RESIDUO
(`self.plDDT_values[residue_number - 1]`) -- con ~7 atomos/residuo en un
PDB de AlphaFold, el "plDDT" del residuo N termina siendo en realidad el
B-factor del residuo ~N/7. Medido en 370 sitios reales: correlacion contra
el B-factor CA real (correcto) = 0.031; con el fix = 0.927 (`avg_plDDT` y
`sasa` estan bien, r=0.9995/0.98 -- el promedio-por-atomo casualmente da
igual).

Recalibrando solo esos 3 tipos con el fix aplicado (mismas muestras):
citrullination 0.657->**0.778** (iguala el paper exacto), s_nitrosylation
0.683->**0.770** (cierra ~53% del gap), glutathionylation (control, ya
buena) 0.873->0.886 sin danarse. **No es techo del modelo** -- aunque los
valores publicados de estos 4 tipos (0.769-0.846) ya son los mas bajos de
los 17, con o sin el fix.

Dos defectos adicionales encontrados que SI midieron pero NO importan en
la practica (documentados, no priorizados): DSSP nunca se ejecuta en
`train_PTM/` pese a que `predict.py:293-295` llama `pose.secstruct()`
(produccion queda ciega a estructura secundaria, 100% "loop"); y
`data_loader.py:148` vs `predict.py:308-316` tienen el orden de columnas
E/H/L invertido (H y E intercambiados). Ablacion real: ninguno de los dos
cambia el AUROC mas de 0.01 -- toda la senal recuperable esta en el fix de
plDDT.

**Coordinacion con el wireado del item 1 (2026-08-01, commit `ae327cd`)**:
los umbrales ya wireados en `Settings.py` se midieron CON el bug de plDDT
activo. Con el fix, citrullination pasa de 0.392 a 0.358 y s_nitrosylation
de 0.501 a 0.486 -- cambios moderados, no catastroficos, pero
`DEEPPTMPRED_CALIBRATED_THRESHOLDS` deberia considerarse **provisional**
para estos 2 tipos (y probablemente para el resto en menor medida, no
medido en los 17) hasta una recalibracion completa con el fix aplicado.

**Salvedades explicitas del subagente**: el fix de plDDT se valido en 3
tipos de 17, no en los 17 -- el resto son extrapolacion razonable, no
medicion directa. `crotonylation` no puede mejorar con mas muestreo
aunque se quisiera: la corrida ya uso TODOS los negativos que existen en
el dataset (23), igual que citrullination con sus 59 positivos -- limite
real de los datos fuente, no de `--n-per-class`.

**Ambas decisiones tomadas y ejecutadas 2026-08-01** (Enzo eligio "aplicar
ahora" en las dos): (a) `n_linked_glycosylation` excluido del consenso de
produccion (`CONSENSUS_EXCLUDED_TYPES` en `ptm_annotation.py`, commit
`05cc4d0`); (b) parche de plDDT aplicado en
`_deepptmpred_runner.py::_load_predict_module` (mismo `05cc4d0`, lee el
B-factor CA real desde `self.pose` de PyRosetta en vez del array
`plDDT_values` mal indexado) y recalibracion completa de los 17 tipos
relanzada con el fix (`--n-per-class 75`, calibracion previa movida a
`DeepPTMPred/data/calibration_STALE_preplddtfix_2026-08-01/`, no borrada).

### Resultado real de la recalibracion post-fix de plDDT (2026-08-01)

```
ptm_type                     auroc     suggested_threshold  n_proteins
gamma_carboxyglutamic_acid   0.989689  0.24807824            44
hydroxylation                0.969697  0.35899624            69
o_linked_glycosylation       0.900667  0.26193630            121
malonylation                 0.890667  0.41699925            145
glutathionylation            0.884444  0.46646200            145
lys_methylation              0.877559  0.43064716            143
ubiquitination                0.876622  0.53213940            149
sumoylation                   0.875315  0.37326753            139
acetylation                   0.878378  0.63506210            140
phosphorylation                0.870244  0.24020174            137
arg_methylation                0.860360  0.34068727            139
crotonylation                  0.815072  0.86312497            34
succinylation                  0.811911  0.50403893            144
citrullination                 0.778079  0.36854228            46
glutarylation                  0.772444  0.47470970            103
s_nitrosylation                0.769730  0.51403310            139
n_linked_glycosylation         0.507200  0.99802846            110
```

**Confirma al 100% el diagnostico del subagente**: los 4 tipos antes
"mediocres" ahora son solidos -- crotonylation 0.607->0.815,
citrullination 0.657->0.778, glutarylation 0.673->0.772, s_nitrosylation
0.683->0.770. El fix de plDDT (validado solo en 3/17 tipos antes de esta
corrida) resulto generalizar a los 4, no solo a los 2 medidos entonces.
**16/17 tipos en rango solido-a-excelente (0.77-0.99)**. Unico tipo que
sigue roto: `n_linked_glycosylation` (0.507, practicamente identico al
0.491 pre-fix) -- confirma que su problema es independiente del bug de
plDDT (SMOTE no-op + dataset 85.6% positivo, ver seccion de arriba), la
exclusion del consenso ya aplicada es la decision correcta.

`Settings.DEEPPTMPRED_CALIBRATED_THRESHOLDS` (wireado 2026-08-01,
`ae327cd` -> reemplazado con estos valores finales) y
`tests/test_settings.py`/`tests/test_ptm_annotation.py` actualizados con
los 17 umbrales de esta tabla. 101 tests pasan.

## Proximos pasos reales

1. ~~Decidir el gap de calibracion `fpr` de DeepMVP~~ — RESUELTO 2026-07-28
   noche (ver seccion arriba). `pipeline.py` corre end-to-end de verdad en
   Camino PDB (572 sitios, 98 con consenso real). En una maquina nueva:
   correr `conda run -n deepmvp python scripts/generate_deepmvp_calibration.py`
   una vez despues de instalar los pesos de DeepMVP, antes del primer uso
   real del pipeline.
2. ~~Repetir/promediar el relax en `ddg_estimate.py`~~ — RESUELTO 2026-07-28
   noche: `--nstruct` (default 3), ver seccion de Fase A / Extension 3
   arriba.
3. ~~Fase A clase 2 (glicosilacion)~~ — IMPLEMENTADA 2026-07-28 noche (ver
   seccion arriba). Default de nucleo biosintetico documentado, no una
   prediccion real de glicoforma. ~~**Mejora: consultar GlyGen antes de caer
   al default generico**~~ — IMPLEMENTADA 2026-07-29: nuevo
   `src/structural/glygen_client.py` (endpoint real descubierto leyendo
   `api.glygen.org/swagger.json`, `POST /protein/detail/{accession}/`,
   verificado con una consulta real en vivo contra P10636/Tau -- 14 sitios
   reales devueltos, confirma el esquema `type`/`start_pos`/`site_category`/
   `glytoucan_ac` usado por el cliente). Wireado como corroboracion PURAMENTE
   INFORMATIVA en `pyrosetta_glycan_patch.py` (`check_glygen_evidence`,
   flag opcional `--uniprot-accession`): no cambia el glicano adjuntado (el
   `glytoucan_ac` que reporta GlyGen identifica un glicano especifico que
   requeriria traducirlo a un nucleo IUPAC construible por PyRosetta -- un
   problema de mapeo real no resuelto por esta mejora, documentado como
   limite explicito), solo informa si el sitio ya tiene evidencia
   experimental conocida. Fallos de red degradan a un aviso, nunca abortan
   el flujo principal. 8 tests nuevos (`tests/test_glygen_client.py`, red
   mockeada -- a diferencia del resto de `src/structural/`, este cliente NO
   requiere `pyrosetta`, solo `urllib`/`json` de la stdlib, asi que si corre
   en la suite principal).
4. ~~Fase A clase 3 (ubiquitinacion/sumoilacion)~~ — IMPLEMENTADA 2026-08-01,
   ver seccion "Fase A clase 3 ... CORRECCION 2026-08-01" arriba. La
   conclusion previa (exige compilar Rosetta completo) era incorrecta:
   `UBQ_GTPaseMover` (via PyRosetta, cero compilacion nueva) cubre ambos
   tipos con el mismo mecanismo, verificado con corridas reales en
   `src/structural/ubiquitin_sumo.py`.
5. Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion —
   deliberadamente sin integrar por ahora (decision 2026-07-26).
6. ~~Calibracion real de DeepPTMPred~~ — IMPLEMENTADA 2026-07-30/31, ver
   seccion "Calibracion real de DeepPTMPred + bug de distribucion
   phi/psi" abajo. ~~**Pendiente decision de Enzo**: si conectar los
   umbrales calibrados a `Settings.py`/produccion~~ — WIREADO 2026-08-01,
   **actualizado el mismo dia** tras encontrar y corregir el bug de plDDT
   (ver seccion "investigacion de n_linked_glycosylation y los 4 tipos
   mediocres" + "Resultado real de la recalibracion post-fix de plDDT"
   abajo): `Settings.DEEPPTMPRED_CALIBRATED_THRESHOLDS` (dict, 17 tipos,
   valores finales post-fix) + `Settings.deepptmpred_threshold_for(tipo)`
   (fallback a `DEEPPTMPRED_MIN_PROBABILITY`=0.5 si el tipo no esta
   calibrado). `ptm_annotation.py` usa el umbral por tipo en los 3 puntos
   donde antes usaba el umbral fijo 0.5 (ambos loops de
   `annotate_pdb_path`); `n_linked_glycosylation` ademas excluido del
   consenso (`CONSENSUS_EXCLUDED_TYPES`) por tener AUROC 0.51 pese al
   umbral bien calibrado -- problema de poder discriminativo del modelo,
   no del corte. 4 tests nuevos (`tests/test_settings.py` +1 en
   `tests/test_ptm_annotation.py` para la exclusion de consenso).

## Auditoria de robustez pre-checkpoint (2026-07-28, noche) -- pulida 2026-07-29

Revision sistematica de todo el codigo (no solo lo ya rastreado arriba)
buscando fragilidad real de cara a un checkpoint estable. Nada de esto
bloqueaba el uso actual del pipeline -- eran mejoras reales de robustez para
"pulir", ordenadas de mayor a menor severidad. **Pulido 2026-07-29: items
1/2/3/4/5/6/7/8 resueltos.** Queda abierto, a proposito, solo el item 10
(licencia de DeepPTMPred) -- requiere contactar a Carlos, no es algo
resoluble desde el codigo. Ver detalle abajo.

**Severidad media (arreglar primero):**

1. ~~**Riesgo real de path traversal, sin sanitizar**~~ — RESUELTO 2026-07-29.
   `DeepPTMPredEngine.run()` construye la carpeta de resultados con
   `base_output_dir / record.accession` -- un join de path CRUDO. `record.accession`
   venia de `structure_parser.py:188` (`accession = path.stem`, el nombre del
   archivo PDB/mmCIF de entrada) y NUNCA se saneaba, a diferencia del accession
   de un registro FASTA (`fasta_parser.py:170-178`, que SI reemplaza `/` y `\`
   explicitamente). Un archivo de entrada llamado literalmente `...pdb` (3
   puntos) produce `accession=".."` (`Path("...pdb").stem == ".."`, verificado
   empiricamente -- la nota original decia `..pdb` con 2 puntos, que en
   realidad da `stem == "."`, no `".."`; corregido aqui), lo que escapa el
   directorio de salida. Arreglo aplicado: nuevo helper
   `structure_parser._sanitize_accession` (mismo criterio que `fasta_parser.py`
   para `/`/`\`, mas fallback a `"UNKNOWN"` para accession vacio/`.`/`..`),
   llamado al derivar `accession` en `parse_structure`. 2 tests nuevos
   (`test_accession_dotdot_se_sanea_evita_path_traversal`,
   `test_sanitize_accession_reemplaza_separadores_de_ruta`).

2. ~~**Cache de features ESM de DeepPTMPred puede servir datos obsoletos sin
   avisar**~~ — RESUELTO 2026-07-29. La clave de cache era solo `protein_id`
   (`{accession}_full_esm.npz`), NO un hash/huella de la secuencia real --
   una secuencia DISTINTA bajo el mismo accession/nombre de archivo (p.ej.
   una version actualizada del PDB con el mismo nombre) reutilizaba en
   silencio el embedding ESM viejo. Arreglo aplicado: nuevo helper
   `_deepptmpred_runner._esm_cache_path` incluye un hash corto (sha256,
   12 hex) de la secuencia real en el nombre del archivo de cache
   (`{protein_id}_{hash}_full_esm.npz`) -- una secuencia distinta ya no
   puede colisionar con una cache existente. 3 tests nuevos en
   `tests/test_deepptmpred_runner.py`.

3. ~~**`BaseEngine.run()` no coincide con las implementaciones reales**~~ —
   RESUELTO 2026-07-29. La interfaz abstracta declaraba
   `run(self, items) -> List[TOut]`, pero `DeepMVPEngine.run` y
   `DeepPTMPredEngine.run` en realidad requerian un segundo parametro
   `output_dir: Path = None` no declarado en el contrato. Arreglo aplicado:
   `output_dir: Optional[Path] = None` anadido a la firma abstracta de
   `src/engines/base_engine.py`. 3 tests nuevos en `tests/test_base_engine.py`.

**Severidad baja (higiene, no urgente):**

4. ~~`requirements.txt` desactualizado~~ — RESUELTO 2026-07-29: comentario
   reescrito, ya no dice "no estan construidos todavia" ni usa la
   numeracion vieja "Fase 3".
5. ~~`ModelLoadError` sin usar~~ — RESUELTO 2026-07-29: clase muerta
   eliminada de `src/utils/exceptions.py` (confirmado por grep que no se
   usaba en ningun otro sitio antes de borrarla).
6. ~~`generate_deepmvp_calibration.py` sobreescribe en silencio~~ —
   RESUELTO 2026-07-29: ahora imprime un aviso explicito con el
   `testing_suffix` usado cuando `site_prediction.tsv` ya existe y va a
   sobreescribirse.
7. ~~El logger de consola solo muestra `WARNING` o mas grave~~ —
   CONFIRMADO INTENCIONAL 2026-07-29 (Enzo). Comportamiento deseado
   (consola limpia); el progreso INFO real sigue disponible en
   `logs/ptm_pipeline.log`. No se cambia codigo. `src/utils/logger_config.py:36`.
8. ~~Sin tests dedicados para `scripts/generate_deepmvp_calibration.py`,
   `src/engines/base_engine.py`, `src/utils/logger_config.py`~~ — RESUELTO
   2026-07-29: `tests/test_generate_deepmvp_calibration.py` (10, cubre
   `trim_window` y `download_all_data` mockeado -- `generate_for_type` sigue
   sin test directo, requiere TensorFlow/`lib` de DeepMVP real, conda env
   `deepmvp`), `tests/test_base_engine.py` (3), `tests/test_logger_config.py`
   (4). `src/structural/*.py` sigue sin tests por diseno (requiere
   `pyrosetta`, ver nota ya existente arriba).

**Salvedad cientifica -- RESUELTA 2026-07-29 (leido el texto real del Methods,
no solo el resumen):**

9. ~~`testing_70` podria ser circular/optimista~~ -- CONFIRMADO NO CIRCULAR.
   Texto real de Methods (DeepMVP, Nature Methods 2025, via PMC12446062,
   texto completo libre):

   > "Ninety per cent of the data was used for training (81%) and
   > validation (9%), and the remaining 10% was used for independent
   > testing."

   > "To control for sequence similarity, peptides in the testing set
   > were filtered to remove those with identities above predefined
   > thresholds (90%, 80% or 70%) compared with peptides in the training
   > and validation sets. This filtering was performed using the
   > clustering tool CD-HIT."

   Es decir: el split primario ES un 10% genuinamente separado del 90%
   train+validation (los pesos shipeados en `models.tar.gz` se entrenaron
   sobre ese 90%, nunca sobre el 10% de test). Los sufijos `70/80/90` NO
   son splits alternativos del dataset -- son el MISMO 10% de test,
   filtrado ADEMAS por similitud de secuencia via CD-HIT contra
   train+validation (70/80/90 = umbral de identidad maximo permitido:
   `testing_70` es el subconjunto MAS estricto, con CUALQUIER peptido de
   ≥70% identidad respecto a train/val ya excluido -- el mas conservador
   de los tres, confirma que `--testing-suffix 70` (default actual del
   script) era ya la eleccion correcta por intuicion, ahora verificada).
   El propio paper reporta que el AUROC se mantiene estable bajo este
   filtro mas estricto ("AUROC values remained stable across all PTMs,
   with changes < 0.03 compared to the original models"), consistente con
   el AUROC 0.90-0.99 medido localmente. **Conclusion: el AUROC medido por
   `generate_deepmvp_calibration.py` es real, no circular** -- seguro de
   citar para TFG/publicacion.

**Sin resolver, ya conocido, re-listado aqui por completitud:**

10. ~~DeepPTMPred no declara licencia~~ -- RESUELTO 2026-07-29: confirmado
    CC BY-NC 4.0 (mismos terminos que el paper) directamente por Junwen
    Wang, autor de correspondencia. Ver seccion de instalacion arriba para
    la cita completa de su respuesta.

## MeToken: corroboracion informativa de TIPO (IMPLEMENTADA 2026-08-01)

Investigado a fondo el 2026-08-01 (subagente Opus) que
`github.com/A4Bio/MeToken` (ICLR 2025) es el motor estructural mas potente
evaluado hasta ahora para PTM -- consume coordenadas backbone reales (N/CA/
C/O) via grafo 3D-kNN + marcos locales por cuaternion, mucho mas rico que
los 4 escalares (SASA/phi/psi/plDDT) que usa DeepPTMPred -- pero el
checkpoint publicado es un CLASIFICADOR DE TIPO en sitios YA CONOCIDOS, no
un detector de sitio (confirmado en `model_interface.py:40`, la clase
"no-PTM" queda enmascarada en entrenamiento cuando `with_null_ptm=0`, como
viene el checkpoint publicado -- y verificado empiricamente contra
`AF-P10636-F1-model_v4.pdb`: predice tipos con alta confianza tambien en
posiciones sin PTM real, prolinas/glicinas). Por eso NO se implemento como
segundo motor de consenso (Decision 2 sigue pausada) sino como corroboracion
PURAMENTE INFORMATIVA del tipo en sitios que el consenso YA acepto
(`pasa_umbral=true`), mismo patron no-decisorio que GlyGen.

### Instalacion real -- verificada, no asumida

- Conda env dedicado `metoken` (`/home/enzo/miniconda3/envs/metoken`, Python
  3.10, CPU-only -- confirmado que esta maquina no tiene `nvidia-smi`).
  `torch==2.13.0+cpu` (indice oficial CPU de PyTorch) + `numpy`/`scipy`/
  `biopython`/`transformers`/`omegaconf`/`tqdm`/`pandas`/`huggingface-hub`/
  `h5py` instalados sin error.
- **`torch_scatter` sin wheel prebuilt para esta combinacion** (confirmado:
  `data.pyg.org/whl/torch-2.13.0+cpu.html` solo lista variantes
  macOS/Windows/aarch64 de `pyg_lib`, ninguna de `torch_scatter` para
  `linux_x86_64` en esa version) -- `pip install --no-build-isolation
  torch_scatter` compilo desde fuente sin error, mas rapido de lo estimado
  de antemano (pocos minutos reales en esta maquina, no los ~15 min
  previstos).
- Repo clonado (`git clone https://github.com/A4Bio/MeToken`, gitignorado,
  mismo patron que `DeepMVP/`/`DeepPTMPred/`). Pesos reales descargados y
  verificados vivos: `https://github.com/A4Bio/MeToken/releases/download/1.0/pretrained_model.zip`
  respondio HTTP 200 real (`content-length: 88260083`, confirmado con
  `curl -sIL` antes de descargar), contiene `checkpoint.ckpt` (36MB, el
  state_dict plano que usa `inference.py` -- confirmado con `torch.load` +
  `model.load_state_dict()` -> `"<All keys matched successfully>"`) y
  `lightning_checkpoint.ckpt` (60MB, checkpoint completo de
  pytorch-lightning con estado de optimizador, sin usar aqui).
- Tokenizer `facebook/esm2_t33_650M_UR50D` (solo el tokenizer, NO el encoder
  ESM-2 completo -- MeToken usa su propio `nn.Embedding` entrenado desde
  cero, `wo_esm`, pese al nombre del atributo) ya estaba cacheado en
  `~/.cache/huggingface/hub/` de una sesion anterior de este mismo dia
  (StackGlyEmbed, proyecto hermano) -- confirmado con `find`, cero descarga
  nueva necesaria para esta corrida. `HF_HUB_OFFLINE=1`/
  `TRANSFORMERS_OFFLINE=1` fijados en `_metoken_runner.py` (mismo patron que
  `stackglyembed_predict_local.py` del proyecto 1).

### Dos bugs reales confirmados corriendo el repo (no asumidos de antemano)

1. **`inference.py:61` llama a `PDB.Polypeptide.three_to_one`**, eliminado
   de Biopython en la version >=1.80. Confirmado real:
   `hasattr(PDB.Polypeptide, 'three_to_one')` -> `False` en Biopython 1.87
   (la version que instala `pip install biopython` hoy). Parcheado en
   `_metoken_runner.py::_patch_three_to_one` -- monkeypatch sobre el modulo
   `Bio.PDB.Polypeptide` ya importado, usando
   `protein_letters_3to1`/`protein_letters_3to1_extended` (los diccionarios
   que SI existen en Biopython moderno) en su lugar. No se edita
   `inference.py` (vendored).
2. **`src/metoken_model.py:213` tiene `device='cuda'` hardcodeado**
   (`codebook_mask = torch.ones(len(codebook), dtype=torch.int32,
   device='cuda')`, dentro de `MeToken_Model.__init__`). Confirmado real
   corriendo SIN parche: revienta con `AssertionError: Torch not compiled
   with CUDA enabled`. Es la UNICA linea de todo `src/` con `device='cuda'`
   hardcodeado (confirmado por grep -- el resto de tensores usan
   `device=x.device`/`device=index.device`, siguiendo el tensor de entrada).
   Parcheado en `_metoken_runner.py::_force_cpu_ones` -- context manager que
   envuelve `torch.ones` SOLO durante la construccion del modelo,
   redirigiendo `device='cuda'`->`'cpu'` cuando CUDA no esta disponible. No
   se edita `src/metoken_model.py` (vendored).

Ambos parches verificados de forma diferencial: corriendo
`examples/Q16613.pdb` (el ejemplo del propio repo, posicion 31) SIN el
parche 2 revienta con el `AssertionError` de arriba; CON ambos parches
reproduce EXACTO el resultado documentado en el propio
`quick_inference.ipynb` del repo (`"PTM type at the position 31 is
Phosphorylation"`) -- confirma que los parches no alteran el comportamiento
numerico del modelo, solo lo hacen correr en CPU.

### Vendorizado

- `src/engines/_metoken_runner.py`: runner standalone (patron
  `_deepptmpred_runner.py`), corre en el venv dedicado `metoken`, recibe un
  PDB de una sola cadena + lista de posiciones 1-based ya aceptadas por el
  consenso, devuelve el tipo con mayor probabilidad de las 24 clases reales
  (excluye indice 0 = "Not a PTM type", la clase null enmascarada en
  entrenamiento, e indice 25 = "in Rare PTM Types", un cubo de tipos raros
  agrupados no interpretable como tipo especifico -- `src/constant.py::
  PTMtype_list` tiene 26 entradas en total, verificado leyendo el archivo
  real).
- `src/engines/metoken_engine.py`: `get_type_corroboration(pdb_path,
  positions, chain_id="A") -> dict[int, dict]`, invoca el runner via
  subprocess. Degradacion NO fatal en todos los modos de fallo (repo no
  instalado, checkpoint ausente, exit code != 0, timeout, `python_bin`
  inexistente, CSV de salida ausente/malformado) -- siempre devuelve `{}` y
  registra un aviso, nunca lanza. 10 tests (`tests/test_metoken_engine.py`,
  `subprocess.run` mockeado, mismo patron que
  `tests/test_deepptmpred_engine.py`).
- `src/engines/ptm_annotation.py::annotate_pdb_path` acepta ahora
  `pdb_path`/`chain_id` opcionales (`None`/`"A"` por defecto -- llamarla sin
  ellos es identico al comportamiento de antes de esta mejora, confirmado
  por los tests ya existentes sin modificar). Si se proveen y
  `Settings.METOKEN_ENABLED`, anade 3 columnas SOLO en filas con
  `pasa_umbral=true`: `metoken_type`, `metoken_probability`,
  `metoken_type_coincide` (`True`/`False` si `tipo_ptm` tiene equivalente
  mapeado en `CANONICAL_TO_METOKEN_TYPE` y coincide o no con
  `metoken_type`, `None` si no hay equivalente conocido). NUNCA toca
  `pasa_umbral`/`consenso` -- doble seguro: el propio `metoken_engine`
  degrada a `{}` sin lanzar, y ademas la llamada esta envuelta en un
  `try/except` en `annotate_pdb_path` por si un fallo no anticipado
  ocurriera en el wiring mismo. 7 tests nuevos en
  `tests/test_ptm_annotation.py` (columnas ausentes sin `pdb_path`, columnas
  solo en filas elegibles, `coincide=True`/`False`/`None`, degradacion a
  `None` si MeToken devuelve `{}`, respeta `METOKEN_ENABLED=False`, no
  tumba la anotacion si el wiring lanza una excepcion inesperada).
- **Hallazgo real durante el testing del wiring, no solo teorico**: `tipo_ptm`
  en las filas DeepMVP-solo (sin match de DeepPTMPred) llega con el nombre
  CRUDO de DeepMVP (`acetylation_k`, `methylation_r`, etc.), no el nombre
  canonico -- si `_add_metoken_corroboration` buscara directamente en
  `CANONICAL_TO_METOKEN_TYPE` sin normalizar primero via
  `DEEPMVP_TO_CANONICAL_TYPE`, `metoken_type_coincide` habria quedado
  `None` incorrectamente para TODA fila DeepMVP-solo (la mayoria de filas
  del reporte). Corregido antes de dar el wiring por terminado -- 2 tests
  (`test_annotate_pdb_path_con_pdb_path_agrega_columnas_solo_en_filas_pasa_umbral`,
  `test_annotate_pdb_path_metoken_type_coincide_false_si_discrepa`) fallaron
  primero SIN la normalizacion, confirmando el bug antes de arreglarlo.
- `pipeline.py::run_fase3_pdb_annotation` pasa `record.chain_pdb_path`
  (el PDB de UNA sola cadena derivado en Fase 1.5, NO `record.pdb_path`
  original que puede tener mas de una cadena) + `record.chain_id` --
  necesario para que las posiciones 1-based que usa MeToken coincidan
  exactamente con `record.sequence`, mismo criterio que ya usa
  `_deepptmpred_runner.py`. Confirmado leyendo `structure_parser.py` que el
  PDB de una sola cadena conserva el nombre de cadena original (nunca lo
  renombra).
- `Settings.py`: `METOKEN_ENABLED` (default `True`, degradacion silenciosa
  si no esta instalado -- no hace falta desactivarlo a mano en una maquina
  sin MeToken), `METOKEN_HOME`, `METOKEN_RUNNER_SCRIPT`, `METOKEN_PYTHON_BIN`,
  `METOKEN_CHECKPOINT`, `METOKEN_TIMEOUT_SECONDS` -- mismo patron que
  `DEEPPTMPRED_*`/`DEEPMVP_*`.

### Verificacion real end-to-end contra Tau (`AF-P10636-F1-model_v4.pdb`)

`pipeline.py` completo (Camino PDB) corrido de verdad TRES veces en esta
sesion con los 3 venvs reales (`deepmvp`, `deepptmpred`, `metoken`)
encadenados via `DEEPMVP_PYTHON_BIN`/`DEEPPTMPRED_PYTHON_BIN`/
`METOKEN_PYTHON_BIN`; los numeros de esta seccion son de la corrida FINAL
(la que ya incluye la normalizacion `phosphorylation_y`->`Phosphorylation`
en `CANONICAL_TO_METOKEN_TYPE`, ver hallazgo abajo). MeToken corrio UNA SOLA
VEZ sobre las 353 posiciones unicas que el consenso acepto
(`pasa_umbral=true`), no una vez por sitio -- ~6s de inferencia real en CPU
para las 353 posiciones (confirmado en el log: `Fase 2 completa` a
`17:56:45`, `Fase 3 completa` a `17:56:51`, incluye la invocacion completa
de MeToken entre medias).

Reporte final (`fasta_outputs/AF-P10636-F1-model_v4_ptm_sites.csv`, 749
filas que pasan el umbral, 116 con consenso real DeepMVP+DeepPTMPred):
`metoken_type` poblado en las 749/749 filas (MeToken predijo algo para
TODAS las posiciones pedidas, ninguna cayo fuera de rango). Desglose real
por `tipo_ptm` (columna `metoken_type_coincide`, solo cuenta filas con
equivalente mapeado en `CANONICAL_TO_METOKEN_TYPE`; `phosphorylation` y
`phosphorylation_y` agregados en una sola fila porque MeToken no distingue
residuo para fosforilacion):

```
tipo_ptm                     coincide/evaluadas   tasa
phosphorylation (+_y)             115/115         100.0%
arg_methylation                    21/26           80.8%
ubiquitination                     37/64           57.8%
acetylation                        26/61           42.6%
sumoylation                         1/59            1.7%
hydroxylation                       0/70            0.0%
gamma_carboxyglutamic_acid          0/57            0.0%
malonylation                        0/53            0.0%
lys_methylation                     0/45            0.0%
succinylation                       0/23            0.0%
o_linked_glycosylation              0/112           0.0%
glycosylation_n                     0/3             0.0%
glutathionylation                   0/3             0.0%
```
(sin equivalente en MeToken, `coincide=None` para 58 filas: crotonylation
50, glutarylation 7, citrullination 1 -- ver `CANONICAL_TO_METOKEN_TYPE`,
estos 3 tipos no tienen clase correspondiente entre las 24 reales de
`src/constant.py::PTMtype_list`.) Conteo total de la columna
`metoken_type_coincide`: `True`=200, `False`=491, `None`=58 (200+491+58=749).
Tasa global sobre las 691 filas evaluables: 200/691 (28.9%).

**Lectura real de estos numeros, no solo la tabla**: MeToken predijo SOLO 6
etiquetas distintas de las 24 posibles en todo Tau (`Phosphorylation` 294,
`Ubiquitination` 272, `Acetylation` 153, `Methylation` 21, `Sumoylation` 7,
`S-palmitoylation` 2) -- fuerte sesgo hacia los 3 tipos de PTM mas comunes
(fosforilacion/ubiquitinacion/acetilacion probablemente dominantes en su
dataset de entrenamiento), casi nunca predice tipos mas raros aunque se le
pregunte directamente en un sitio que otro motor SI acepto para ese tipo:

- **Fosforilacion: acuerdo perfecto (115/115, incluye los 6 sitios
  `phosphorylation_y` -- MeToken no distingue S/T de Y, una sola clase
  "Phosphorylation" para las 3, verificado en `src/constant.py`)**. Unico
  tipo donde MeToken corrobora consistentemente.
- **Sumoilacion: casi nunca corrobora (1/59)** -- pero el patron es
  biologicamente coherente, no ruido aleatorio: de los 59 sitios de
  sumoilacion (residuo K, confirmado), MeToken predice `Ubiquitination` en
  37 y `Acetylation` en 21 (las otras 2 modificaciones de Lys que si estan
  entre sus 6 etiquetas usadas, ver distribucion completa abajo). Ubiquitina
  y SUMO comparten exactamente el mismo mecanismo quimico (enlace
  isopeptidico Lys-Gly, ver `src/structural/ubiquitin_sumo.py` de Fase A
  clase 3) -- una confusion microambiente-nivel entre las modificaciones mas
  parecidas quimicamente del set (todas requieren K), no un fallo
  arbitrario. Ejemplo real: posicion 657 (K657, `score_deepmvp=0.9997`,
  `score_deepptmpred=0.79`, consenso=True) es el UNICO sitio de sumoilacion
  donde MeToken si dice "Sumoylation" (probabilidad 0.426); en el resto
  (44, 67, 87, 130, 189, 204...) dice "Ubiquitination" o "Acetylation" con
  probabilidades 0.44-0.62.
- **Hidroxilacion: 0/70, con un hallazgo mas fuerte que "no corrobora"**:
  las 70 posiciones son TODAS prolina (`residuo_wt='P'`, confirmado),
  correcto para hidroxilacion -- pero MeToken predice `Phosphorylation` (49
  veces) o `Ubiquitination` (21 veces), tipos que biologicamente exigen
  S/T/Y o K, NUNCA P. El top-1 de MeToken en estos 70 sitios ignora la
  identidad quimica del propio residuo -- mas que "MeToken no conoce
  hidroxilacion", sugiere que el microambiente de estas prolinas en Tau
  (proteina intrinsecamente desordenada, plDDT medio 49.34 en este modelo
  AlphaFold, ver Extension 3 arriba) no produce una senal fuerte y el
  modelo cae de vuelta a sus clases mayoritarias.
- **`S-palmitoylation` (tipo sin ningun equivalente en el pipeline, ninguno
  de los 17 tipos de DeepMVP/DeepPTMPred lo cubre) aparecio 2 veces**, ambas
  en sitios que el pipeline etiqueta `glutathionylation` (posiciones 608 y
  639) -- otra confusion quimicamente plausible (ambas son modificaciones de
  cisteina).
- **Acetilacion (42.6%) y ubiquitinacion (57.8%, contando SOLO cuando
  `tipo_ptm='ubiquitination'`, no las 37 confusiones de sumoilacion de
  arriba) son parcialmente corroboradas** -- ambas comparten LYS como
  residuo objetivo con `arg_methylation`/`lys_methylation`/`sumoylation`,
  consistente con que MeToken distinga razonablemente bien "algo le pasa a
  esta Lys" pero no siempre acierte cual de las 4 modificaciones especificas
  de Lys es la correcta.

**Conclusion practica**: MeToken corrobora con alta fiabilidad SOLO
fosforilacion en este caso real; para el resto de tipos, sus discrepancias
son mayoritariamente explicables por confusion entre modificaciones
quimicamente relacionadas (mismo residuo objetivo) en vez de ruido -- util
como corroboracion informativa tal como esta diseñado (nunca decide
`pasa_umbral`/`consenso`), pero NO deberia leerse como "MeToken desmiente
estos sitios": su cobertura de tipo esta fuertemente sesgada hacia las PTMs
mas comunes de su propio dataset de entrenamiento, un limite real del
checkpoint publicado (no del wiring de este proyecto).

118 tests (`pytest tests/`, subio de 101 -- 10 de
`tests/test_metoken_engine.py` + 7 de `tests/test_ptm_annotation.py`; el
resto de `src/engines/`/`src/utils/` no cambiaron). `MeToken/` gitignorado
(mismo patron que `DeepMVP/`/`DeepPTMPred/`/`MTPrompt-PTM/`).
