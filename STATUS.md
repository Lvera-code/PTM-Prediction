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
decision escrita en `01-Proyectos/PTM-Prediction/Decisiones/2026-08-01-fase-a-clase3-ubiquitinacion-sumoilacion-implementada.md`.

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
5. ~~Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion~~
   — IMPLEMENTADA 2026-08-01, ver seccion "StackGlyEmbed: corroboracion
   informativa de N-GLICOSILACION" abajo. Reusa el venv/pickles ya
   instalados en `B-Cell-Epitope-Prediction` como recurso externo (sin
   importar codigo entre proyectos). Corroboracion puramente informativa,
   nunca decide `pasa_umbral`/`consenso` -- verificada con una corrida real
   contra los 3 candidatos de N-glicosilacion de Tau.
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
   **REVERTIDO 2026-08-01** (Enzo pidio explicitamente limpiar esto):
   nivel del handler de consola subido a INFO -- el progreso (ejecucion
   de motores, routing por accession) ya es visible en vivo, no solo en
   el archivo rotativo. Test actualizado
   (`test_console_handler_nivel_info`, antes `_warning`). Commit `fef3e3e`.
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

Reporte final (`outputs/AF-P10636-F1-model_v4_ptm_sites.csv`, 749
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

## StackGlyEmbed: corroboracion informativa de N-GLICOSILACION (IMPLEMENTADA 2026-08-01)

Item 5 de "Proximos pasos" (deliberadamente sin integrar desde el
2026-07-26) ejecutado: cross-validacion de `n_linked_glycosylation` con
StackGlyEmbed del proyecto 1. Motivo (ver seccion "investigacion de
n_linked_glycosylation" arriba): DeepPTMPred no tiene NINGUN poder
discriminativo real para este tipo (AUC 0.495 en las metricas de
entrenamiento del propio repo, no arreglable reentrenando -- ya excluido del
consenso via `CONSENSUS_EXCLUDED_TYPES`). StackGlyEmbed
(`github.com/GaryChan-lab/StackGlyEmbed`) es un tercer motor
INDEPENDIENTE de arquitectura (ProteinBERT + ESM-2 650M + ProtT5 apilados ->
meta-clasificador SVM), especializado UNICAMENTE en N-glicosilacion --
mismo patron no-decisorio que MeToken/GlyGen (nunca cambia
`pasa_umbral`/`consenso`, solo corrobora un candidato ya propuesto).

### Recurso externo reusado -- verificado real, no asumido

A diferencia de MeToken/DeepPTMPred (repos clonados DENTRO de este
proyecto), StackGlyEmbed ya estaba instalado y funcionando de verdad en el
proyecto HERMANO `B-Cell-Epitope-Prediction` (proyecto 1, independiente por
decision 2026-07-26 -- nunca se importa codigo entre proyectos, pero SI se
reusan recursos externos pesados ya instalados). Verificado con `ls` real
antes de asumir nada:

- Venv `B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed`
  (Python 3.10.20, symlink `python -> python3.10`) -- existe real.
- `B-Cell-Epitope-Prediction/StackGlyEmbed/prediction/` -- 3
  `power_transformer_*.sav` + `base_layer_pickle_files/` (30 modelos base:
  10 SVM + 10 XGB + 10 KNN, mas `SVM_meta_layer.sav`, ~480MB total) --
  existen reales, tamanos de archivo consistentes con pesos reales
  (10-33MB cada uno).
- ProtT5 (`Rostlab/prot_t5_xl_half_uniref50-enc`, reusado de
  `scipion-chem-tmbed/tmbed_src/tmbed/models/t5/`) -- `model.safetensors`
  de 3.2GB presente, real.

### Vendorizado (NO se importa codigo del proyecto hermano, decision 2026-07-26)

- `src/engines/_stackglyembed_runner.py`: adaptacion PROPIA (no una copia
  importada) de la logica real ya verificada en
  `B-Cell-Epitope-Prediction/src/engines/stackglyembed_predict_local.py`
  (extraccion de features ProteinBERT+ESM-2+ProtT5 + prediccion via el
  stack de clasificadores ya entrenado). Simplificado respecto al script
  del hermano: como aqui siempre se evalua UNA proteina por invocacion
  (mismo patron que `_deepptmpred_runner.py`/`_metoken_runner.py`), el CLI
  recibe `--sequence`/`--positions` directamente (sin el formato
  intermedio `dataset.txt` multi-proteina del repo original). Mismo
  criterio de chunking para proteinas completas (ESM-2 por bloques de
  1024aa, ProtT5 por 8797aa) y misma ventana de 15 residuos para el
  promedio de ESM-2 -- verificado que la posicion 1-based del secuon
  coincide exactamente con la convencion de `sequence`/`posicion` del resto
  del pipeline (secuencia completa, sin offsets de ventana que convertir).
- `src/engines/stackglyembed_engine.py`: `get_nglyco_corroboration(sequence,
  positions, result_dir=None, filename_prefix="") -> dict[int, dict]`,
  invoca el runner via subprocess sobre
  `Settings.STACKGLYEMBED_PYTHON_BIN` (venv del proyecto hermano).
  Degradacion NO fatal en todos los modos de fallo (venv ausente, runner
  ausente, pickles ausentes, exit code != 0, timeout, CSV de salida
  ausente/malformado) -- siempre devuelve `{}` y registra un aviso, nunca
  lanza. 11 tests (`tests/test_stackglyembed_engine.py`, `subprocess.run`
  mockeado, mismo patron que `tests/test_metoken_engine.py`).
- `src/engines/ptm_annotation.py`: **a diferencia de MeToken (requiere un
  PDB, Camino PDB unicamente), StackGlyEmbed solo necesita `sequence` --
  ya un parametro EXISTENTE en ambas funciones (`annotate_fasta_path`/
  `annotate_pdb_path`), asi que aplica a AMBOS caminos sin necesitar ningun
  parametro nuevo relacionado con PDB.** Ambas funciones aceptan ahora
  `enable_stackglyembed: bool = False` (backward-compatible: llamarlas sin
  el parametro es identico al comportamiento de antes). Si es `True` y
  `Settings.STACKGLYEMBED_ENABLED`, para cada fila con `tipo_ptm` en
  {`n_linked_glycosylation`, `glycosylation_n`} (cubre tanto el nombre
  canonico de DeepPTMPred como el crudo de DeepMVP -- nunca se fusionan en
  consenso, ver `CONSENSUS_EXCLUDED_TYPES`) y `pasa_umbral=True`, anade 3
  columnas informativas: `stackglyembed_veredicto`
  (`'Glicosilado'`/`'No glicosilado'`), `stackglyembed_score` (probabilidad
  cruda del meta-SVM) y `stackglyembed_coincide` (`True` si el veredicto es
  `'Glicosilado'` -- corrobora el candidato, `False` si lo contradice; sin
  la ambiguedad de "tipo sin equivalente" de MeToken, porque StackGlyEmbed
  solo predice N-glicosilacion). NUNCA toca `pasa_umbral`/`consenso` --
  mismo doble seguro que MeToken (el propio engine degrada a `{}` sin
  lanzar, mas un `try/except` en el wiring). 8 tests nuevos en
  `tests/test_ptm_annotation.py`.
- `pipeline.py`: `run_fase3_fasta_annotation`/`run_fase3_pdb_annotation`
  pasan `enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED` a
  `annotate_fasta_path`/`annotate_pdb_path` -- unico cambio en el
  orquestador, sin nuevos parametros de PDB/cadena (no hacen falta).
- `Settings.py`: `STACKGLYEMBED_ENABLED` (default `True`, degradacion
  silenciosa si el recurso externo no esta disponible),
  `STACKGLYEMBED_PYTHON_BIN` (default apunta al venv REAL del proyecto
  hermano, verificado con `ls`), `STACKGLYEMBED_RUNNER_SCRIPT` (el propio
  vendorizado de este proyecto), `STACKGLYEMBED_MODELS_DIR` (default apunta
  a `B-Cell-Epitope-Prediction/StackGlyEmbed/prediction`, verificado),
  `STACKGLYEMBED_T5_MODEL_PATH`, `STACKGLYEMBED_ESM_MODEL_NAME`,
  `STACKGLYEMBED_TIMEOUT_SECONDS` (900, mismo valor que usa el proyecto
  hermano para la misma carga en frio de 3 modelos).

### Verificacion real end-to-end contra Tau (`AF-P10636-F1-model_v4.pdb`)

`get_nglyco_corroboration` corrido de verdad (no mockeado) contra el venv
REAL del proyecto hermano, sobre la secuencia completa de Tau (758
residuos, la misma ya extraida por Fase 1.5 en una corrida previa de este
proyecto) y los 3 candidatos REALES de `n_linked_glycosylation` que el
pipeline ya habia aceptado (`pasa_umbral=true`) en la corrida completa mas
reciente (`outputs/AF-P10636-F1-model_v4_ptm_sites.csv`): posiciones
484 (secuon `NAT`), 676 (`NIT`), 727 (`NVS`). Carga en frio real de los 3
embedders (ProteinBERT + ESM-2 650M + ProtT5) + inferencia: **52.5s** de
tiempo real (`time`, no estimado).

Resultado real (`stackglyembed.csv` persistido):

```
posicion  sequon  veredicto        score
484       NAT     No glicosilado   0.2110
676       NIT     No glicosilado   0.3943
727       NVS     No glicosilado   0.1217
```

Tambien verificado el WIRING completo (no solo el engine aislado): llamando
`annotate_pdb_path` real con estos 3 candidatos como `deepmvp_df` y
`enable_stackglyembed=True` reproduce EXACTOS los mismos 3 scores
(`stackglyembed_coincide=False` en las 3 filas) -- confirma que la
integracion en `ptm_annotation.py` invoca el motor real correctamente, no
solo en aislamiento.

**Lectura real de estos numeros**: StackGlyEmbed NO corrobora ninguno de
los 3 candidatos que DeepMVP propuso para Tau (los 3 scores caen del lado
"No glicosilado", aunque 676 con 0.39 queda mas cerca del limite que los
otros dos). Esto es informativo, no una prueba definitiva en ningun
sentido: por un lado, refuerza la cautela ya documentada sobre
`n_linked_glycosylation` en este proyecto (el propio DeepPTMPred no tiene
poder discriminativo aqui); por otro, Tau es una proteina intrinsecamente
desordenada (plDDT medio 49.34, ver Extension 3 arriba) y no hay evidencia
experimental externa (GlyGen no devolvio N-glicosilacion confirmada para
P10636 en las consultas ya hechas, ver seccion de Fase A clase 2) que
confirme cual de los tres motores acierta en este caso concreto -- ningun
motor tiene la ultima palabra, StackGlyEmbed es una corroboracion adicional
puramente informativa, exactamente como se diseño.

136 tests (`pytest tests/`, subio de 118 -- 11 de
`tests/test_stackglyembed_engine.py` + 8 de `tests/test_ptm_annotation.py`;
el resto de `src/engines/`/`src/utils/` no cambiaron). El recurso externo
(`B-Cell-Epitope-Prediction/StackGlyEmbed/`) vive fuera de este repo, nada
nuevo que gitignorar aqui.

## Decision 2 (segundo motor real de consenso, Camino PDB) -- CERRADA 2026-08-06 (EMNGly)

No confundir con [[MeToken]] (implementado arriba como corroboracion
informativa de TIPO, no como motor de consenso -- limitacion real
encontrada: el checkpoint publicado es un clasificador de tipo en sitios
ya conocidos, no un detector de sitio).

### CoNglyPred (candidato original, 2026-08-01) -- confirmado DEFINITIVAMENTE muerto

El candidato original para el rol de segundo motor SIMETRICO en consenso
(como DeepMVP+DeepPTMPred, pero acotado a `n_linked_glycosylation`, el
unico tipo con motor muerto) era **CoNglyPred**
(`github.com/whm242446/CoNglyPred`, Proteomics 2025) -- graph transformer
sobre PDBs de AlphaFold2 + DSSP + co-atencion con ESM-2. Verificado
exhaustivamente 2026-08-01 que no tenia pesos publicados en ningun sitio.
Enzo envio un correo el 2026-08-01 a Shaoping Shi (autora de
correspondencia, Universidad de Nanchang) pidiendo los pesos -- **sin
respuesta a los 5 dias**. Re-verificado 2026-08-06 via API de GitHub
directamente (no por busqueda): 0 releases/tags/issues, ultimo push
2024-08-15, CERO archivos `.pth`/`.pt`/`.pkl` en todo el arbol pese a que su
propio README instruye cargar `best.pth`. Confirmado genuinamente muerto,
no solo lento en responder -- descartado.

### EMNGly (reemplazo adoptado, 2026-08-06)

Investigacion de reemplazo delegada a un subagente Opus (permiso permanente,
ver vault) mas re-verificacion propia clonando el repo real. **EMNGly**
(`github.com/StellaHxy/EMNgly`, Hou et al., Bioinformatics 39(11):btad650,
2023) es el unico candidato con pesos REALES y descargables que preserva la
propiedad de diseno ya decidida al descartar MTPrompt-PTM como reemplazo de
DeepPTMPred (el segundo motor debe usar estructura 3D real, no solo
secuencia):

- Arquitectura: ESM-1b (`site_emb`+`local_emb`, 1280+1280, ver
  `model/get_esm_embedding.py`) + **MIF** (Masked Inverse Folding de
  Microsoft, `model/MIF/`, embedding estructural REAL de 256-dim sobre
  coordenadas de backbone N/CA/C del PDB, licencia BSD-2 permisiva --
  verificada en el repo oficial de Microsoft, la copia vendorizada en
  EMNgly perdio el LICENSE) -> SVM (`sklearn.svm.SVC`, 2816 features).
- Pesos verificados A NIVEL DE BYTES 2026-08-06 (no de README): MIF
  (`mif.pt`, 13.8MB) bundled en el repo; SVM (`N-GlyDE.pickle`, 36MB) en una
  carpeta publica de Google Drive, descargado con una peticion HTTP Range
  real sin autenticacion -- confirmado pickle protocolo 4,
  `sklearn.svm._classes.SVC`, `n_features_in_=2816` (coincide exacto con
  1280+1280+256), `_sklearn_version=1.1.1`.
- Se eligio `N-GlyDE.pickle` (no `N-GlyAltas_classifier.pkl`, tambien
  disponible en la misma carpeta) porque esta entrenado sobre el benchmark
  con NEGATIVOS RESTRINGIDOS AL SEQUON N-X-[S/T] -- el regimen correcto para
  este rol: un modelo entrenado contra negativos sin restringir (como usan
  MIND-S/PDeepPP, candidatos evaluados y descartados por este motivo) puede
  acertar solo aprendiendo el motivo que este pipeline YA aplica via
  `_nglyco_window`, sin aportar una segunda opinion real.
- Benchmark publicado (N-GlyDE test set): EMNGly MCC 0.736 vs LMNglyPred
  0.717, DeepNGlyPred 0.605, N-GlyDE 0.499, NetNGlyc 0.265.
- Sin LICENSE declarado en ningun repo (el paper SI es CC BY 4.0 real,
  confirmado 2026-08-07 via el XML de PMC/PMC10627407, no cubre el codigo)
  -- a diferencia de CoNglyPred, esto NO es bloqueante: los pesos ya estan
  descargables. Correo REDACTADO 2026-08-07 (envio programado 2026-08-10,
  lunes) a los autores de correspondencia reales -- Yaojun Wang
  (`wangyaojun@cau.edu.cn`, China Agricultural University) y Shiwei Sun
  (`dwsun@ict.ac.cn`, ICT-CAS), ambos verificados via el mismo XML de PMC
  (la pagina de Oxford Academic esta detras de un challenge de Cloudflare,
  bloquea fetch directo) -- mismo patron que funciono con DeepPTMPred, no
  bloquea produccion. Ver vault, addendum 2026-08-07 en
  `2026-08-06-decision2-cerrada-emngly-reemplaza-conglypred.md`, para el
  texto completo del correo.

**Implementado** (ver `src/config/settings.py` bloque `EMNGLY_*`,
`src/engines/_emngly_runner.py`, `src/engines/emngly_engine.py`,
`src/engines/ptm_annotation.py::_apply_nglyco_consensus`): consenso REAL
de 3 motores (DeepMVP+EMNGly+StackGlyEmbed) para
`n_linked_glycosylation`/`glycosylation_n` en Camino PDB --
`pasa_umbral` = al menos 1 de los motores disponibles pasa (generaliza la
regla OR del resto de tipos), `consenso` = al menos
`Settings.NGLYCO_CONSENSUS_MIN_ENGINES` (default 2) pasan. StackGlyEmbed
se PROMUEVE de corroboracion informativa (rol que mantiene tal cual en
Camino FASTA) a motor de consenso real en Camino PDB. Fix proactivo de
alineamiento incluido en el runner: `structure_emb` de MIF esta indexado
por NUMERO DE RESIDUO DEL PDB (no por orden ATMSEQ) -- el runner traduce
via la tabla `position_mapping` de Fase 1.5 (`fasta_position`/`pdb_seqid`)
en vez de asumir que ambas numeraciones coinciden (cierto solo para PDBs de
AlphaFold2 sin huecos, como los que usa el propio dataset de entrenamiento
de EMNGly -- ver docstring de `_emngly_runner.py` para el analisis completo,
incluye tambien la aclaracion de que la convencion de indices BOS-inclusive
de `site_emb`/`local_emb` NO es un bug pese a la apariencia inicial).

**2026-08-07: los 2 go/no-go checks PASARON, `Settings.EMNGLY_MIN_PROBABILITY=0.5`
ya NO es provisional.** Pesos reales descargados en esta maquina por primera
vez (ESM-1b 7.8GB + companero de regresion de contactos, SVM `N-GlyDE.pickle`
35.9MB via HTTP Range de Google Drive, confirmado `sklearn.svm.SVC` real con
`n_features_in_=2816`) y corridos de punta a punta.

**4 bugs reales encontrados y arreglados en el camino** (ninguno tocando
codigo vendorizado, mismo criterio que el resto del proyecto):
1. `_emngly_runner.py` solo anadia `EMNgly/model` a `sys.path` -- correcto
   para el propio `from MIF.sequence_models... import`, pero
   `MIF/sequence_models/*.py` (`pdb_utils.py`, `pretrained.py`, etc.) hacen
   imports BARE de `sequence_models.xxx` (no `MIF.sequence_models.xxx`),
   necesitan `EMNgly/model/MIF` en `sys.path` TAMBIEN -- confirmado corriendo
   el import real, no solo leyendo codigo. El script original
   `get_mif_embedding.py` nunca expuso este problema porque ademas hace
   `sys.path.append("./MIF")` Y depende de que su propio directorio ya este
   en `sys.path` (comportamiento automatico de un script, no de un modulo
   importado via subprocess).
2. Esa misma cadena de imports depende de `wget` (paquete pip, dependencia
   transitiva de `trRosetta_utils.py` nunca usada por este proyecto) --
   ausente, no documentado en el `environment.yml` de EMNgly.
3. `pip install scikit-learn==1.1.1` sin pinnear numpy resuelve numpy 2.x
   por defecto -- binariamente incompatible con el wheel de sklearn 1.1.1
   (`ValueError: numpy.dtype size changed`). Fix: pinnear `numpy==1.23.5`
   despues.
4. `fair-esm` (paquete pip, no vendorizado) llama `torch.load(path,
   map_location="cpu")` sin `weights_only=False` -- desde PyTorch 2.6 el
   default cambio a `True` y el checkpoint de fair-esm no pasa el allowlist
   estricto (`argparse.Namespace`). Parche acotado en
   `_ESMEmbeddingExtractor.__init__` (unico choke point que carga el
   checkpoint), mismo patron que el monkeypatch de
   `_deepptmpred_runner.py::_load_predict_module`.

**Check 1 (alineamiento) PASO**: Alpha-1-antitrypsin (P01009) tiene 3 sitios
N-glicosilacion reales confirmados en GlyGen (`reported_with_glycan`,
consulta real via `glygen_client.py`) en numeracion UniProt 70/107/271 --
numeracion PDB (mature protein, offset -24) 46/83/247. PDB real usado: 1QLP
("2.0 Angstrom structure of intact alpha-1-antitrypsin", RCSB), con 22
residuos N-terminales genuinamente ausentes (`REMARK 465`, confirmado) --
el caso de "PDB con huecos" que le importa a este fix. Fase 1.5 real
(`parse_structure`) sobre 1QLP: `fasta_position` 24/61/225 -> `pdb_seqid`
46/83/247, los 3 confirmados `residue_letter='N'`. Verificado a 2 niveles:
(a) `MIF.sequence_models.pdb_utils.parse_PDB`'s `wt[pdb_seqid-1]` = `'N'`
para los 3 -- Y se demostro que el indexado NAIVE sin la traduccion
(`wt[fasta_position-1]`, el bug que este fix previene) da `'F'`/`'G'` en 2 de
los 3 casos reales, no `'N'` -- bug real, no hipotetico. (b) El runner
completo (ESM-1b + MIF + SVM reales) sobre esos 3 sitios: probabilidades
0.887/0.949/0.729 (los 3 por encima del umbral 0.5). Control negativo real:
2 secuones rotos (`N-x-P`, biologicamente nunca glicosilados) en la misma
proteina, posiciones 82/345 -> probabilidades 0.526/0.286 (separacion real,
aunque 82 queda justo en el borde del umbral -- esperado, un modelo real no
es perfecto).

**Check 2 (MCC) PASO Y SUPERO el numero publicado**: `MCC=0.8197` (AUC
0.9631, especificidad 0.914, sensibilidad 0.906) sobre 301 sitios reales
evaluables del set independiente de N-GlyDE (`NGLYDE_independent.txt`, via
`dukkakc/DeepNGlyPred` -- el propio repo de EMNgly solo documenta la fuente,
no la empaqueta), vs. el `MCC=0.736` publicado en el paper. 146/447 filas
(33%) quedaron fuera: 4/86 proteinas sin modelo AlphaFold DB (78 filas) mas
fallos de alineamiento puntuales en las 82 restantes. Estructura real
usada: AlphaFold DB, API `https://www.alphafold.ebi.ac.uk/api/prediction/`
(NUNCA construir la URL `AF-{acc}-F1-model_v4.pdb` a mano -- confirmado real
2026-08-07 que AlphaFold DB ya sirve v6 y v4 da 404 en el 100% de los casos
probados, mismo hallazgo y mismo fix que
`scripts/generate_deepptmpred_calibration.py` ya aplico antes). Superar el
numero publicado con margen (no solo "acercarse") es razonable, no
sospechoso: los modelos AlphaFold v6 usados aqui son mas recientes/precisos
que los disponibles cuando el paper de EMNGly se publico (2023) -- una
mejora legitima de la fuente estructural, no fuga de datos (el SVM ya
estaba entrenado, congelado, nunca se reentreno con este set).

`predict.py::get_scores` tiene un bug real confirmado (no arreglado en el
vendorizado, ver `scripts/verify_emngly_nglyde_mcc.py::get_scores` para el
puerto corregido): llama `get_scores(label_y, predict_y[0])` -- `predict_y[0]`
es un escalar (la probabilidad de la PRIMERA fila), no la lista completa.

Scripts nuevos (permanentes, reproducibles, gitignored solo los datos que
descargan): `scripts/prepare_emngly_nglyde_structures.py` (Fase 1.5 real por
proteina, venv `cnb_pipeline`), `scripts/verify_emngly_nglyde_mcc.py`
(embeddings + SVM + scores, venv `.venv-emngly`).

Instalacion (venv dedicado `emngly`, pin `scikit-learn==1.1.1` -- ver
hallazgo de version arriba, README.md tiene el detalle completo incluyendo
los 4 bugs reales de instalacion documentados arriba):
```bash
git clone https://github.com/StellaHxy/EMNgly
python3 -m venv .venv-emngly
.venv-emngly/bin/pip install fair-esm torch "scikit-learn==1.1.1" scipy pandas numpy tqdm wget
.venv-emngly/bin/pip install "numpy==1.23.5"
# ESM-1b (~7.4GB) + companero de regresion de contactos, descarga manual (ver README.md)
# SVM: https://drive.usercontent.google.com/download?id=1hbnEtHHXTGnQAFm-cCHMj3pWQiAYAUsw&export=download&confirm=t
```

## Fase A conectada al pipeline principal (2026-08-03)

Revierte la decision 2026-07-27 de que D (`apply_workflow_filter`) no ruta a
Extension 3/Fase A porque esas fases no existian todavia -- motivo: demo
completa de punta a punta a Carlos Oscar Sorzano el 2026-08-10, se pidio
explicitamente conectar Fase A (no solo el nucleo de anotacion) para esa
fecha. Decision completa en el vault,
`01-Proyectos/PTM-Prediction/Decisiones/` (buscar la nota del 2026-08-03).

**Problema real identificado antes de escribir codigo**: modelar TODOS los
sitios que pasan el umbral es computacionalmente inviable -- un caso real
como Tau acepta ~572 sitios, y cada modelado estructural (ddG con nstruct=3,
conjugacion, refinado de glicano) tarda minutos, no segundos. Decision
explicita de Enzo: **top-N por tipo** (default 1, `Settings.FASE_A_TOP_N_PER_TYPE`)
en vez de todos los sitios o solo bajo demanda manual -- representativo y
demostrable en una corrida real automatica.

### Piezas nuevas

- `src/structural/fase_a_dispatch.py`: enruta un sitio (`pdb_path`,
  `position`, `ptm_type`) al modulo real correspondiente segun 3 clases
  mutuamente excluyentes (ver docstring del modulo para el detalle completo):
  clase 1 (5 tipos, `pyrosetta_ptm_patch.py` + `ddg_estimate.py`), clase 2 (2
  tipos, `pyrosetta_glycan_patch.py`), clase 3 (2 tipos, `ubiquitin_sumo.py`).
  Los otros 8/17 tipos devuelven `estado="sin_soporte_fase_a"` sin tocar
  PyRosetta. Cada llamada asume un proceso PyRosetta fresco -- las 3 clases
  inicializan PyRosetta con flags incompatibles entre si (confirmado leyendo
  los 3 modulos: `pyrosetta.init()` solo aplica las flags de su PRIMERA
  llamada por proceso), por eso nunca se mezclan clases en el mismo proceso
  (ver `_fase_a_runner.py`, un sitio por subprocess).
- Arreglo previo necesario: `ddg_estimate.py` importaba `pyrosetta_ptm_patch`
  con un import BARE (`from pyrosetta_ptm_patch import ...`), que solo
  funcionaba invocando el script como standalone desde su propio directorio
  -- cambiado a `from src.structural.pyrosetta_ptm_patch import ...` (mismo
  patron que `pyrosetta_glycan_patch.py` ya usaba para `glygen_client`), para
  que `fase_a_dispatch.py` pueda importar los 3 modulos de forma fiable.
- `src/engines/_fase_a_runner.py`: runner standalone (subprocess, un sitio
  por proceso, mismo patron que `_deepptmpred_runner.py`) -- inserta la raiz
  del repo en `sys.path` el mismo (a diferencia de los runners de motores
  externos, este importa codigo PROPIO del proyecto, no un repo vendorizado,
  asi que necesita el insert explicito). Serializa el resultado a JSON.
- `src/engines/fase_a_engine.py`: `FaseAEngine(BaseEngine)` + dataclass
  `FaseASiteRequest`. Subprocess por sitio via `Settings.FASE_A_PYTHON_BIN`
  (reusa el MISMO conda env `deepptmpred` ya instalado para DeepPTMPred --
  mismo PyRosetta real, no un env nuevo). Degradacion NO fatal en todos los
  modos de fallo (interprete/runner ausentes, exit code != 0, timeout, JSON
  ausente/malformado) -- mismo patron que `stackglyembed_engine.py`/
  `metoken_engine.py`: un sitio que falla nunca tumba el resto del barrido.
- `src/engines/ptm_annotation.py::select_fase_a_candidates`: selecciona el
  top-N por tipo (prioriza `score_deepptmpred`, fallback `score_deepmvp`)
  restringido a `Settings.FASE_A_SUPPORTED_PTM_TYPES` (constante de datos
  pura, sin pyrosetta -- `fase_a_dispatch.py` valida en tiempo de import que
  coincide exactamente con la union real de los 3 `SUPPORTED_PTM_TYPES`,
  falla alto si alguna vez divergen).
- `pipeline.py`: `run_fase3_pdb_annotation` ahora devuelve tambien el
  DataFrame en memoria (no solo `report_path`); nuevo paso
  `run_fase_a_pdb_modeling` despues de Fase 3 en el Camino PDB -- reescribe
  el CSV final con las columnas `fase_a_estado`/`fase_a_clase`/`fase_a_ddg`/
  `fase_a_glycan_tree`/`fase_a_glygen_evidencia`/`fase_a_conjugation_metrics`/
  `fase_a_output_pdb` para TODAS las filas aceptadas (no solo las
  seleccionadas): las no seleccionadas quedan `fase_a_estado="no_seleccionado"`,
  documentado explicitamente en vez de una columna vacia ambigua. Camino
  FASTA no cambia (Fase A requiere PDB).

18 tests nuevos (154 total): `tests/test_fase_a_engine.py` (11, subprocess
mockeado, mismo patron que StackGlyEmbed/MeToken), 7 en
`tests/test_ptm_annotation.py` para `select_fase_a_candidates` (logica pura
de pandas, sin pyrosetta, tests reales no mockeados), mas un test de
integracion en `tests/test_pipeline_fase1.py` mockeando `FaseAEngine.run`
para verificar el wiring completo del reporte final.

### Corrida real end-to-end (2026-08-03): 7/8 candidatos modelados, 1 hallazgo real nuevo

Verificado con una corrida real completa (Camino PDB: Fase 1.5->2->3->A) sobre
`DeepPTMPred/data/AF-P10636-F1-model_v4.pdb` (Tau), con `FASE_A_PYTHON_BIN`/
`DEEPPTMPRED_PYTHON_BIN` apuntando al mismo conda env `deepptmpred` real (no
mockeado en ningun punto). 749 sitios totales en el reporte final, 8
candidatos seleccionados por `select_fase_a_candidates` (no 9: ningun sitio
de `n_linked_glycosylation` supero su umbral calibrado de 0.998 -- el tipo se
salta limpiamente, confirmado, sin candidato que modelar para el).

```
tipo_ptm                     posicion  fase_a_estado  detalle
acetylation                  465       modelado       ddG = +50.72 (ref2015_cart)
lys_methylation               460       modelado       ddG = -16.90
phosphorylation                529       modelado       ddG = +42.52
gamma_carboxyglutamic_acid     319       modelado       ddG = +124.50
hydroxylation                  499       error          ver hallazgo abajo
o_linked_glycosylation         498       modelado       glicano = core_1_O-glycan
sumoylation                    215       modelado       conjugacion real (5 metricas de torsion)
ubiquitination                  755       modelado       conjugacion real (5 metricas de torsion)
```

Tiempos reales: ~7-8min por sitio de clase 1 (ddG, nstruct=3 + relax final
para el PDB de salida), ~20-40s para clase 3 (conjugacion), ~1-2min para
clase 2 (glicano). Corrida completa (incluyendo Fase 1/2/3 previas: DeepMVP,
17 invocaciones de DeepPTMPred, MeToken, StackGlyEmbed) en ~43 minutos.

**Hallazgo real nuevo: `hydroxylation` falla siempre via
`add_variant_type_to_pose_residue`, no es un bug de este proyecto.**
Investigado en profundidad, no solo registrado: el patch real de Rosetta
para esto (`pro_hydroxylated_case1.txt`/`case2.txt`,
`TYPES HYDROXYLATION1`/`HYDROXYLATION2`) SI esta cargado (confirmado via
introspeccion de `ResidueTypeSet.patches()`, 122 patches totales, ambos
presentes) y su `BEGIN_SELECTOR` (`PROPERTY PROTEIN`, `NAME3 PRO`,
`HAS_ATOMS 2HG`) SI se cumple en la pose real (el residuo PRO499 tiene el
atomo `2HG` explicito, confirmado listando sus atomos). Pese a esto,
`add_variant_type_to_pose_residue` con `VariantType.HYDROXYLATION`,
`.HYDROXYLATION1` y `.HYDROXYLATION2` (las 3 variantes existen como enum
real) fallan las 3 con el mismo error:

```
ERROR: Unable to find desired residue 'PRO' with variant 'HYDROXYLATION1'.
Attempted to add target variant(s) to ResidueType using both ResidueType
base name 'PRO' and base ResidueType.
```

Verificado que esto NO es especifico de este PDB/pose (mismo resultado
probando limpio, sin ningun otro patch aplicado antes). La causa raiz exacta
(por que `add_variant_type_to_pose_residue` no resuelve una combinacion
base+variante que si esta declarada y cuyo selector se cumple) no se
investigo mas a fondo por presupuesto de tiempo -- **no se aplico ningun
workaround sin poder validar su correccion quimica**, mismo criterio ya
establecido para el patch roto de ubiquitina (2026-07-28/08-01): mejor
documentar un limite real y confirmado que forzar una solucion no verificada.
`FaseAEngine`/`fase_a_dispatch.py` degradan esto correctamente
(`estado="error"`, mensaje real capturado, el resto del barrido continua sin
verse afectado) -- el wiring en si funciono exactamente como se diseño; el
limite es de Rosetta/`pyrosetta_ptm_patch.py`, no del wiring nuevo. Pendiente
para una sesion futura con mas tiempo: los otros 4 tipos de clase 1
(acetylation, lys_methylation, phosphorylation, gamma_carboxyglutamic_acid)
funcionaron sin problema -- este limite es especifico de hidroxilacion.

### Robustez: ddG std ya no se descarta (2026-08-03, mismo dia)

Auditoria post-demo-prep encontro que `ddg_estimate.estimate_ddg` ya
calculaba la desviacion estandar entre las `nstruct` trayectorias
independientes de relax (WT y parcheado), pero `fase_a_dispatch._run_class1`
solo propagaba el minimo (`wt_score`/`mut_score`), descartando `wt_scores`/
`mut_scores` antes de que la incertidumbre llegara al reporte final -- un
ddG se mostraba mas seguro de lo que realmente era. Arreglado: nuevo
`Settings.FASE_A_RESULT_TEMPLATE` (unica fuente de verdad de las claves de
un resultado de Fase A, compartida entre `fase_a_engine.py`/
`fase_a_dispatch.py`/`_fase_a_runner.py`, evita que las 4 rutas de salida
queden desalineadas) incluye ahora `ddg_std`/`wt_score_std`/`mut_score_std`;
`pipeline.py` anade `fase_a_ddg_std` al CSV final. 154 tests siguen
pasando (actualizado `_EMPTY_EXTRA` en `test_fase_a_engine.py`).

### Robustez, puntos 5/6 (2026-08-04): modo batch + tests dedicados de `src/structural/`

Continuacion del plan de robustez/produccion post-demo-prep (puntos 1-4:
LICENSE, ddG std, `requirements.txt` fijado, CI basico, ver arriba).

**Punto 5 (modo batch)**: `pipeline.py --input` acepta ahora un directorio
ademas de un archivo unico. `src/utils/input_router.py` ya tenia
`route_inputs` (plural) preparado desde su diseño original pero nunca se
habia conectado a la CLI -- `main()` solo llamaba a `route_input` (singular)
sobre un unico archivo. Extraido `run_single_input()` (misma logica de
`main()` de antes, factorizada para que el modo de un solo archivo y el modo
batch reusen exactamente el mismo camino, sin duplicarlo) y añadido
`_run_batch()`: descubre archivos con extension reconocida en el directorio
(no recursivo, mismo alcance plano que `inputs/`), corre cada uno,
**un archivo que falla no detiene el resto del batch** (mismo criterio de
degradacion no fatal que `FaseAEngine`/StackGlyEmbed/MeToken aplican
por-sitio), y escribe `batch_summary.csv` en `output_dir` (columnas
`archivo`/`estado`/`detalle`). Codigo de salida 0 solo si TODOS los archivos
completaron sin error; un directorio sin ningun archivo reconocido tambien
falla explicitamente (1), no en silencio. 3 tests nuevos
(`tests/test_pipeline_batch.py`): batch multi-archivo OK, degradacion
parcial (1 OK + 1 invalido, ambos se procesan, exit code 1), directorio
vacio de inputs reconocibles.

**Punto 6 (tests dedicados de `src/structural/`)**: ninguno de los 5 modulos
de Fase A (`fase_a_dispatch.py`, `ddg_estimate.py`, `pyrosetta_ptm_patch.py`,
`pyrosetta_glycan_patch.py`, `ubiquitin_sumo.py`) tenia tests propios --
solo se ejercitaban indirectamente via `test_fase_a_engine.py`/
`test_pipeline_fase1.py`, que mockean `FaseAEngine`/`FaseASiteRequest`
entero sin llegar nunca a la logica interna de estos modulos. Los 5 modulos
solo importan `pyrosetta` dentro de funciones (nunca a nivel de modulo, por
diseño -- ver docstring de `src/structural/__init__.py`), lo que los hace
importables sin PyRosetta instalado; donde una funcion hace un import
incondicional de `pyrosetta.*` ANTES de su propia logica de validacion
(`apply_ptm_patch`/`attach_glycan`/`_run_class2`/`_run_class3`/
`estimate_ddg`), se stubea un arbol minimo de modulos en `sys.modules` (via
`monkeypatch.setitem`) en vez de mockear a nivel de subprocess -- mismo
criterio de "mockear en el limite mas bajo posible" que el resto del
proyecto. 36 tests nuevos, ninguno requiere PyRosetta/conda envs:

- `test_fase_a_dispatch.py` (7): routing a clase 1/2/3 segun `ptm_type`,
  tipo sin soporte (`estado="sin_soporte_fase_a"`, ningun modulo de PyRosetta
  se inicializa), excepcion real de un submodulo traducida a `estado="error"`
  sin propagar, calculo de `ddg_std` (suma en cuadratura de las desviaciones
  WT/mutante) verificado contra `statistics.pstdev`, corroboracion GlyGen
  opcional (clase 2) solo se consulta si se pasa `uniprot_accession`,
  consistencia `SUPPORTED_PTM_TYPES` vs `Settings.FASE_A_SUPPORTED_PTM_TYPES`.
- `test_ubiquitin_sumo.py` (11): las 4 excepciones reales de `_validate_target`
  (tipo no soportado, pose multicadena, numeracion PDB no secuencial --
  el "segundo hallazgo real" documentado en el modulo --, residuo objetivo
  no es Lisina -- Rosetta mismo no lo valida, `runtime_assert` comentada en
  su .cc real), caso valido por cada tipo soportado, `_lyx_params_path`
  (encuentra/no encuentra el params file relativo a `pyrosetta.__file__`),
  y una regresion real: los 2 PDB de referencia empaquetados
  (`ubiquitin_reference.pdb`/`sumo1_reference.pdb`) existen en disco.
- `test_pyrosetta_glycan_patch.py` (5): las 4 ramas de `check_glygen_evidence`
  (red no disponible, sin sitio reportado, sitio con glicano especifico,
  sitio sin glicano especifico) -- esta funcion no necesita ningun stub de
  pyrosetta, solo mockea `glygen_client.lookup_site` en el mismo limite que
  ya usa `test_glygen_client.py`.
- `test_pyrosetta_ptm_patch.py` (9): las 2 excepciones de `apply_ptm_patch`
  (tipo sin `VariantType` nativo, residuo real no coincide con el esperado
  -- en ambos casos `add_variant_type_to_pose_residue` NO se llega a invocar,
  verificado explicitamente), caso valido para los 5 tipos soportados,
  consistencia `PTM_VARIANT_MAP`/`PTM_TARGET_RESIDUE`/`SUPPORTED_PTM_TYPES`.
- `test_ddg_estimate.py` (4): `_best_of_n_relax` devuelve el minimo de
  `nstruct` trayectorias independientes (una pose fresca por trayectoria, no
  reusada) y solo aplica el parche cuando `ptm_type is not None` (rama WT);
  `estimate_ddg` end-to-end con `nstruct=3` y `nstruct=1`, verificado que usa
  el minimo por estado (no el promedio) para el ddG final.

193 tests en total (154 + 3 batch + 36 de `src/structural/`), todos pasando.

### Punto 7 (2026-08-04): causa raiz real de `hydroxylation` encontrada y arreglada

Cierra el ultimo pendiente de Fase A clase 1 (ver "Fase A conectada al
pipeline principal", 2026-08-03, arriba: "hallazgo real nuevo... no
investigado mas a fondo por presupuesto de tiempo"). Investigado a fondo
via un subagente Opus con PyRosetta real (conda env `deepptmpred`), no
solo lectura de codigo.

**Causa raiz**: `pro_hydroxylated_case1.txt`/`case2.txt` son los UNICOS 2
patches de todo `fa_standard` que declaran `SET_BASE_NAME` (`HYP`/`0AZ`,
confirmado via `grep -rl SET_BASE_NAME` sobre el directorio completo de
patches). Un patch con `SET_BASE_NAME` construye un `ResidueType` BASE
nuevo, no una variante del residuo original -- verificado por introspeccion
en caliente: `HYP`/`0AZ` son `is_base_type()==True` con
`base_name()=="HYP"`/`"0AZ"`, nunca `"PRO"`. `add_variant_type_to_pose_residue`
busca dentro de la familia del `base_name()` ACTUAL del residuo (`"PRO"`) --
como `HYP`/`0AZ` nunca aparecen ahi, la busqueda falla con el error ya
documentado, pese a que el patch esta cargado y su selector se cumple. Es
estructural (propiedad fija del patch), ningun flag de init lo cambia. Los
otros 4 tipos de esta clase no tienen este problema porque sus patches NO
declaran `SET_BASE_NAME`.

**Arreglo, en `src/structural/pyrosetta_ptm_patch.py`**: nuevo
`PTM_BASE_NAME_MAP = {"hydroxylation": "HYP"}`; `apply_ptm_patch` reemplaza
el `ResidueType` completo por `HYP` (case1, 4-hidroxi-L-prolina trans -- el
UNICO producto posible de las prolil-4-hidroxilasas PHD/EGLN sobre
sustratos no colagenicos como HIF-1a, que es el caso real que este
pipeline valida; ni DeepMVP ni DeepPTMPred distinguen case1/case2) via
`replace_pose_residue_copying_existing_coordinates` en vez de
`add_variant_type_to_pose_residue` -- mismo patron ya establecido en
`ubiquitin_sumo.py` de usar la API de Rosetta correcta para una clase de
residuo que la API generica no puede resolver.

**Verificado empiricamente** (no solo argumentado), contra PRO499 de
`AF-P10636-F1-model_v4.pdb` (Tau): la geometria ideal reconstruida desde el
icoor del patch coincide con case1 (CG-OD1 1.424A, OD1-HOD 0.970A, chi4
-127.11 grados -- NO con case2's +90.01 grados, confirma la estereoquimica
trans/4R correcta), backbone sin mover (CA/CG desplazamiento 0.0000A), sin
clashes (contacto intra-residuo mas cercano a OD1: 2.16A),
`relax_neighborhood` corre sin error sobre el resultado parcheado, y sin
regresion en los otros 4 tipos de clase 1 (acetylation/lys_methylation/
phosphorylation/gamma_carboxyglutamic_acid, re-verificados en la misma
sesion). 2 tests nuevos en `tests/test_pyrosetta_ptm_patch.py` (11 en
total en ese archivo): la nueva ruta de reemplazo de ResidueType, y que
`PTM_BASE_NAME_MAP` solo cubre `hydroxylation`. 195 tests en total, todos
pasando.

**Cierra el plan de robustez/produccion de Fase A clase 1 al 100%**: los 5
tipos soportados (phosphorylation, acetylation, hydroxylation,
gamma_carboxyglutamic_acid, lys_methylation) funcionan sin ningun
limite conocido pendiente.

### Punto 8 (2026-08-04): panel de validacion biologica (7 proteinas reales adicionales a Tau)

Hasta ahora el pipeline solo se habia corrido de punta a punta contra UNA
proteina real (Tau/MAPT). Investigado por un subagente Opus (criterio ya
autorizado: verificar una afirmacion cientifica contra fuente primaria
antes de fijarla en produccion) que sourcing real de sitios PTM
documentados en literatura para un panel mas amplio, con la regla dura de
este proyecto de nunca fabricar datos cientificos.

**Metodo real de verificacion** (no memoria/aproximacion): cada sitio del
panel viene de una anotacion de UniProt con evidencia `ECO:0000269`
(experimental) y un PMID confirmado via NCBI eutils -- se descarto
explicitamente toda anotacion `ECO:0000250` ("por similitud", no
experimental en humano). El propio subagente encontro y descarto 2 PMIDs
que recordaba de memoria pero que resultaron ser papers no relacionados al
verificarlos -- no estan citados en el panel final.

**Hallazgo real importante que cambia el diseño**: las URLs `_v4` de
AlphaFold estan muertas (404) -- AlphaFold DB esta en v6, hay que usar
`AF-<accession>-F1-model_v6.pdb`. Verificado tambien que **la numeracion de
AlphaFold coincide exactamente con la numeracion canonica de UniProt**
(confirmado corriendo `src.utils.structure_parser.parse_structure` sobre
los 7 PDBs reales descargados: longitud de cadena = longitud UniProt en
los 7 casos) -- por eso el panel usa exclusivamente AlphaFold, no PDBs
experimentales. **3 trampas de numeracion reales encontradas y evitadas**
(verificadas via registros DBREF de PDBs experimentales reales, NO
usadas en el panel final mas que como advertencia): histonas H3/H4 (PDB
experimental = UniProt - 1, el Met inicial se recorta), protrombina (PDB
experimental = UniProt - 43, numeracion de cadena madura), EPO (PDB
experimental = UniProt - 27) -- documentadas en el docstring de
`biological_panel.py` para que nadie las reintroduzca sin querer al añadir
un PDB experimental en el futuro.

**Panel final, 7 proteinas + 1 control negativo real** (`src/validation/
biological_panel.py`, PDBs reales en `inputs/validation_panel/`,
descargados y verificados -- cada residuo del ground truth se comprueba
contra la secuencia real del PDB en `tests/test_biological_panel.py`, 30
tests, todos pasando):
- **p53/TP53** (P04637, 393 aa): fosforilacion/acetilacion/metilacion de
  Lys y Arg/ubiquitinacion/sumoilacion -- 7 sitios tier A (S15/S20/S46/
  S392 fosfo, K382 acetil, K370/K372 metil, todos multi-PMID) + 14 tier B.
- **Histona H3.1** (P68431, 136 aa) e **Histona H4** (P62805, 103 aa):
  las marcas clasicas de metilacion/acetilacion/fosforilacion (H3K4/K9/
  K27/K36 me, H3K9/K14/K18/K23 ac, H3S10/S28 ph, H4K5/K8/K12/K16 ac,
  H4K20 me, H4S47 ph -- todas tier A, decadas de evidencia independiente)
  mas crotonilacion/succinilacion/glutarilacion/citrulinacion tier B (MS a
  gran escala). Nombres de campo (`H3K4` etc) = UniProt - 1, documentado
  explicitamente para no confundir con la posicion real usada.
- **Protrombina/F2** (P00734, 622 aa): **el sitio mas fuerte del panel** --
  10/10 residuos gamma-carboxiglutamato confirmados independientemente por
  registros MODRES de un cristal real (6C2W) ademas de la secuenciacion
  original, mas 3 sitios de N-glicosilacion. Unica fuente real de
  `gamma_carboxyglutamic_acid` en el panel.
- **HIF-1a** (Q16665, 826 aa): unica fuente real de `hydroxylation` (P402/
  P564/N803, tier A pese a 1-4 PMIDs cada uno por ser papers landmark
  universalmente aceptados) mas sumoilacion/ubiquitinacion. pLDDT bajo
  (60.8, proteina intrinsecamente desordenada fuera del dominio bHLH-PAS)
  -- esperado, no es un fallo de datos. K532 acetilacion DELIBERADAMENTE
  EXCLUIDA (contestada en la literatura, el subagente no pudo verificar
  las citas de refutacion que recordaba).
- **EPO** (P01588, 193 aa): 4 sitios de N/O-glicosilacion, tier A pese a 1
  PMID cada uno por venir de una caracterizacion quimica DIRECTA (mapeo de
  peptidos de la glicoproteina purificada, Lai et al. 1986), no de un
  screening.
- **KIT ligand/SCF** (P21583, 273 aa): **control negativo real** -- Asn-97
  es un sequon N-X-S perfectamente valido que UniProt confirma
  explicitamente como NO glicosilado en ninguna de las 2 isoformas
  conocidas (mismo PMID que documenta los sitios SI glicosilados de la
  misma proteina) -- distingue un motor con especificidad real de uno que
  solo empareja el motivo N-X-[S/T].

Cobertura: **9/9 tipos de Fase A** y 15/17 tipos del nucleo (todos menos
malonylation/glutathionylation, que se dejaron fuera del panel final para
controlar alcance -- ver el reporte completo del subagente para 2 proteinas
opcionales, GAPDH/CLIC1, que las cubririan si se quiere ampliar despues).

**Construido y verificado (no solo diseñado)**:
1. Los 7 PDBs reales descargados de `alphafold.ebi.ac.uk` y verificados
   parseables por `src.utils.structure_parser.parse_structure` (Fase 1.5
   real, no simulada) -- longitud exacta esperada en los 7 casos.
2. `tests/test_biological_panel.py` (30 tests): cada `GroundTruthSite`
   coincide con el residuo real de esa posicion en el PDB descargado (el
   chequeo que habria detectado cualquiera de las 3 trampas de numeracion
   si el panel las hubiera heredado por error), todos los tipos PTM son de
   los 17 soportados, cobertura de los 9 tipos de Fase A confirmada,
   sin duplicados, al menos 1 control negativo real.
3. `scripts/validate_biological_panel.py`: corre Fase 1.5->2->3 REAL
   (sin Fase A -- mide calidad de anotacion/consenso, no modelado
   estructural, que ya se valido por separado) sobre cada proteina del
   panel, compara contra el ground truth, reporta recall tier A/B por
   separado (nunca mezclados) y marca explicitamente cualquier falso
   positivo sobre el control negativo. `tests/test_validate_biological_panel.py`
   (4 tests, motores mockeados, mismo criterio que test_pipeline_fase1.py).

**Primera corrida REAL end-to-end completada (Histona H4, 2026-08-04)**:
motores DeepMVP/DeepPTMPred de verdad, sin mock -- 13 minutos totales
(mayormente DeepPTMPred: 17 tipos, uno por subprocess, ESM-2 domina el
tiempo independientemente del tamaño de la proteina, no es proporcional a
sus 103 residuos). Resultado real: **Tier A 6/6 (100%)**, **Tier B 7/9
(78%)**, sin falsos positivos en el control negativo (H4 no tiene ninguno
propio, ver `kit_ligand_scf`). Los 2 misses de tier B (`S52 phosphorylation`,
`K13 succinylation`) son exactamente el tipo de sitio que el propio panel
documenta como mas debil (screening MS unico, estequiometria real
desconocida) -- consistente con la expectativa, no una señal de alarma.
Efecto secundario real observado, no bloqueante: MeToken fallo con
`ModuleNotFoundError: No module named 'omegaconf'` en el env `cnb_pipeline`
que orquesta el pipeline (MeToken deberia correr en su propio conda env
`metoken` via `METOKEN_PYTHON_BIN`, no configurado en esta corrida) --
degrado exactamente como esta diseñado (corroboracion opcional, log de
aviso, consenso/pasa_umbral sin afectar). Fase A intento modelar 6
candidatos sin `FASE_A_PYTHON_BIN` apuntando al env con PyRosetta (0/6
modelados) -- tampoco afecta el recall medido, que es sobre Fase 3, no
Fase A.

**Corrida de las 6 proteinas restantes completada (2026-08-04)** -- la
corrida en background se interrumpio a mitad de `kit_ligand_scf` (el
control negativo, ultima proteina de la cola) por una caida de la
conexion; 5/6 (p53, histone_h3, prothrombin, hif1a, epo) habian terminado
para entonces con su `_ptm_sites.csv` real ya en disco. Se reanudo solo la
proteina faltante (`scripts/validate_biological_panel.py --only
kit_ligand_scf`, no hizo falta repetir las otras 5) y se calculo el recall
final leyendo los 6 `_ptm_sites.csv` ya generados directamente (misma
logica de `_recall_by_tier`/`_negative_control_report`, sin re-correr el
pipeline sobre las 5 ya completas).

**Resultado real, panel completo (7/7 proteinas)**:

| Proteina | Tier A | Tier B |
|---|---|---|
| p53 | 7/7 (100%) | 14/15 (93%) |
| histone_h3 | 11/11 (100%) | 9/10 (90%) |
| histone_h4 | 6/6 (100%) | 7/9 (78%) |
| prothrombin | 10/13 (77%) | -- |
| hif1a | 5/8 (62%) | 2/4 (50%) |
| epo | 3/4 (75%) | -- |
| kit_ligand_scf | 3/5 (60%) | -- |
| **Global** | **45/54 (83%)** | **32/38 (84%)** |

**Control negativo (Asn-97 de kit_ligand_scf/SCF, sequon N-X-S valido pero
UniProt confirma que NO esta glicosilado): correctamente NO aceptado, 0
falsos positivos.** Distingue un motor con especificidad real de uno que
solo empareja el motivo N-X-[S/T] -- ver diseño del control en
`biological_panel.py`.

Lectura real de los dos casos mas debiles (prothrombin 77%, hif1a
62%/50%): ambas son proteinas grandes con dominios estructuralmente
dificiles (protrombina tiene multiples dominios plegados independientes;
HIF-1a es mayormente intrinsecamente desordenada, pLDDT medio 60.8 fuera
del dominio bHLH-PAS, como ya documentado arriba) -- consistente con la
expectativa de que el pipeline funciona mejor en regiones bien
estructuradas, no una señal de un bug nuevo. No investigado mas a fondo
en esta sesion (analisis por-sitio de los misses queda como trabajo futuro
si se quiere mejorar el recall en proteinas grandes/desordenadas).

229 tests en total, todos pasando.

### Auditoria final pre-demo (2026-08-04): bug real de terminal-nueva encontrado y arreglado

Enzo pidio verificar a fondo que absolutamente todo estuviera listo antes de
la primera sesion con Carlos. Se encontro un bloqueante real y reproducible,
no hipotetico: `DEEPMVP_PYTHON_BIN`/`DEEPPTMPRED_PYTHON_BIN`/
`FASE_A_PYTHON_BIN`/`METOKEN_PYTHON_BIN` caen por defecto a `sys.executable`
(el python que corre `pipeline.py`) si no se exportan -- a diferencia de
`STACKGLYEMBED_PYTHON_BIN` (que ya tiene un default hardcodeado real en
`settings.py`), estos 4 dependian de que Enzo corriera el `export` manual de
README.md en la sesion de terminal activa. Ese export nunca quedo persistido
en `~/.bashrc` (a diferencia de `BEPIPRED_PYTHON_BIN` del proyecto 1, que si
esta ahi). **Confirmado real**: correr el comando exacto del "Quick Start"
del README contra un PDB en una terminal nueva simulada (`bash -ic`, sin
ningun export previo) rompe de inmediato en el primer motor
(`ModuleNotFoundError: No module named 'pyteomics'`, DeepMVP corriendo bajo
el python de `cnb_pipeline`).

**Arreglado**: las 4 lineas `export` se agregaron a `~/.bashrc` (mismo
patron que `BEPIPRED_PYTHON_BIN`), fuera del repo (no es un archivo
versionado del proyecto). **Verificado con una corrida real completa**,
terminal nueva simulada, Camino PDB contra Tau: Fase 1.5 -> Fase 2 (434
DeepMVP + 1003 DeepPTMPred) -> MeToken -> StackGlyEmbed -> Fase 3 (749/1022
pasan umbral, 116 con consenso) -> Fase A (**8/8 candidatos modelados con
exito**, incluyendo `hydroxylation` -- reconfirma en produccion el fix del
Punto 7) -> reporte final generado sin ningun error. Cero excepciones en
todo el log.

Tambien revisado y descartado como bloqueante real: docstring desactualizado
en `deepptmpred_engine.py` (dice "NO PROBADO TODAVIA contra el entorno real",
escrito 2026-07-27, contradicho por docenas de corridas reales desde
entonces) -- cosmetico, no funcional, no arreglado en esta sesion.

**Conclusion de la auditoria**: no queda ningun bloqueante tecnico para las
sesiones con Carlos. El unico item real abierto en todo el proyecto es
Decision 2 (CoNglyPred, segundo motor para `n_linked_glycosylation`) --
esperando respuesta de Shaoping Shi (email 2026-08-01 + seguimiento
2026-08-04, sin respuesta aun en ninguno de los dos hilos).

### Punto 9 (2026-08-04, REVERTIDO el mismo dia): imagen Docker -- descartada

Se construyo un Dockerfile de 4 stages (base/deepmvp/deepptmpred/full) como
continuacion del plan de robustez/produccion post-demo-prep. Causo 2
incidentes reales de salud del host el mismo dia (el build de `full` colgo
la VM de WSL2 por agotamiento de RAM a las 13:07-13:14; por separado, las
imagenes ya construidas -- 27GB -- casi agotaron el disco con solo 11GB
libres en el host). Al revisar la dinamica real del proyecto (Enzo presenta
en su propio ordenador con Carlos al lado, itera en vivo hasta que el lo
declara listo para Scipion -- ver vault), se confirmo que Docker no aporta
nada a ese flujo: los conda envs ya instalados directamente en el host
cubren "vivir en mi ordenador" por completo, y la integracion final a
Scipion se hace via conda envs (igual que el resto de plugins scipion-chem-*
de Enzo), no via contenedores. Decision de Enzo: eliminar `Dockerfile`,
`.dockerignore`, `docker/` y las imagenes ya construidas -- no volver a
proponer containerizar salvo que Carlos o la integracion a Scipion lo
requieran explicitamente. El wheel de PyRosetta (licencia academica) se
conservo, movido a la raiz del repo (gitignorado) para no repetir el
problema de mirrors si se necesita reinstalar en el futuro.

## Analisis de coherencia biologica (2026-08-07): 3 mejoras implementadas, 2 diferidas

A peticion de Enzo, analisis completo del proyecto (robustez de produccion +
coherencia biologica del workflow, equivalente para este proyecto a las
recomendaciones de Carmen Elena Gómez para
[[project_carmen_elena_gomez_feedback_2026-07-30|BCell-Epitope-Prediction]]).
5 recomendaciones concretas encontradas, ninguna un bug de codigo -- huecos
de comunicacion de alcance/interpretacion. Implementadas las 3 de mayor
valor/menor costo; las otras 2 (ambigüedad de tipo de cadena en Fase A de
ubiquitinacion/SUMOilacion, falta de especificidad de quinasa en
fosforilacion) quedan documentadas en README.md "Alcance e interpretacion"
como limitaciones de alcance, sin implementacion planeada.

**Cambio 1 -- disclaimer de interpretacion**: ningun motor modela la via
biosintetica real del sustrato (co-expresion/co-localizacion de la enzima)
-- todos predicen CAPACIDAD de modificarse desde secuencia/estructura, nunca
ocurrencia real en una celula/tejido/condicion. `pipeline.py::INTERPRETATION_DISCLAIMER`
se imprime una vez al final de cada corrida (no se escribe dentro del CSV --
cambiaria su esquema, rompiendo lectores existentes).

**Cambio 2 -- evidencia de via secretora (N-glicosilacion)**: nuevo
`src/structural/uniprot_localization_client.py` (stdlib-only, mismo patron
que `glygen_client.py`), consulta real a `rest.uniprot.org` (endpoint
`/uniprotkb/{accession}.json`, verificado en vivo, no asumido). Columna
`via_secretora_evidencia` (True/False/None), puramente informativa (nunca
decide `pasa_umbral`/`consenso`, mismo patron que MeToken/GlyGen), aplicada
solo a `n_linked_glycosylation`/`glycosylation_n` (deliberadamente NO a
`o_linked_glycosylation` -- dos vias biologicas distintas, O-GlcNAc vs
mucina, que este cliente no distingue). Toggle: `Settings.SECRETORY_PATHWAY_CHECK_ENABLED`
(default True).

3 bugs reales encontrados verificando contra la API real (no asumidos):
1. UniProt devuelve HTTP 400 (no solo 404) para accessions con formato
   invalido -- el caso MAS COMUN aqui, ya que `accession` normalmente viene
   del stem del archivo de entrada, casi nunca un ID UniProt real. Ambos
   codigos se tratan como "sin datos" (`None`), nunca lanzan.
2. La keyword `'Membrane'` sola es demasiado amplia -- descartada tras un
   falso positivo real contra GAPDH/P04406 (proteina citoplasmatica/nuclear
   con reportes de asociacion periferica de membrana).
3. Contar CUALQUIER keyword como "hay datos de localizacion" (sin filtrar
   por `category`) daba un falso `False` en vez del `None` correcto cuando
   UniProt solo tenia keywords no relacionadas con localizacion (p.ej.
   `'3D-structure'`, categoria real `'Technical term'`) -- corregido
   filtrando por `category == 'Cellular component'`.

Verificado en vivo contra 5 casos reales: P01588/EPO (True, secretada),
P04406/GAPDH (False, citoplasmatica/nuclear, localizacion real conocida),
`1qlp`/`NOTREAL999` (None, accession no reconocido).

**Cambio 3 -- aviso de competencia/crosstalk entre PTMs**: varios tipos
modifican el MISMO grupo quimico de un residuo y son mutuamente excluyentes
en una misma molecula/instante (nunca modelado antes de esta mejora, cada
tipo/posicion se puntuaba independiente). Grupos reales, conservadores,
verificados contra `residuo_wt` de cada fila (nunca asumidos solo por el
tipo) en `src/engines/ptm_annotation.py::_PTM_COMPETITION_GROUPS`: acilo-
lisina (acetilacion/ubiquitinacion/sumoilacion/metilacion-Lys/malonilacion/
glutarilacion/succinilacion/crotonilacion), tiol de cisteina (S-nitrosilacion/
glutationilacion), guanidino de arginina (metilacion-Arg/citrulinacion),
hidroxilo de Ser/Thr (fosforilacion/O-glicosilacion, hipotesis "Yin-Yang").
Columna `ptm_crosstalk_aviso`, puramente informativa. Toggle:
`Settings.PTM_CROSSTALK_CHECK_ENABLED` (default True).

**Riesgo real evitado durante el testing**: muchos tests preexistentes de
`test_ptm_annotation.py` (secciones StackGlyEmbed/EMNGly) usan
`glycosylation_n`/`pasa_umbral=True` -- sin mockear, el Cambio 2 habria hecho
llamadas de red REALES a UniProt en cada uno de ellos durante la suite
principal (viola la convencion de este proyecto de nunca golpear una API
real en tests). Arreglado con un fixture `autouse` (`_mock_uniprot_lookup`)
que mockea el limite de red a `None` por defecto para todo el archivo. 23
tests nuevos (14 en `test_ptm_annotation.py`, 9 en el nuevo
`tests/test_uniprot_localization_client.py`), 276 tests en total, todos
pasando, suite completa en ~37s (sin red real).

**Ítem de licencia aclarado en esta sesion (no confundir)**: el email ya
respondido por Junwen Wang ("I confirm that the GitHub code follows the
same CC BY-NC terms") es sobre **DeepPTMPred** (autores Yong Liu/Junwen
Wang), ya cerrado desde 2026-07-29 (ver seccion de licencias arriba). La
licencia de **EMNGly** (autores distintos: Yaojun Wang/Shiwei Sun, ver
seccion de licencias arriba para sus emails reales) sigue siendo el unico
item de compliance real abierto en todo el proyecto -- correo REDACTADO
en esta misma sesion, envio programado para el 2026-08-10 (lunes). No
bloqueante (los pesos de EMNGly ya son descargables sin depender de esa
respuesta).

## Fase A / Fase 3c ELIMINADA del alcance -- integracion a Scipion confirmada (2026-08-10)

Reunion de feedback con Carlos sobre este proyecto (demo del 2026-08-10 ya
realizada). Veredicto: **Fase A/3c (modelado estructural real via PyRosetta
-- ddG, glicosilacion, ubiquitinacion-sumoilacion, todas las secciones
anteriores de este documento fechadas 2026-07-28/2026-08-03) queda
eliminada por completo del alcance del proyecto**, no diferida. Nucleo
(Fase 1 -> 1.5 -> 2 -> 3 -> 3b: DeepMVP+DeepPTMPred, consenso, via
secretora, Kinase Library, MeToken, EMNGly+StackGlyEmbed para
N-glicosilacion, competencia entre PTMs) **confirmado feature-complete** --
"no hay que agregar nada mas". Carlos confirma ademas que **el proyecto se
integrara en Scipion** (mismo patron que
[[project_carlos_scope_decision_2026-07|la integracion del proyecto 1]]).
Escritura completa de la decision, con lo que queda pendiente decidir sobre
orden de trabajo:
`01-Proyectos/PTM-Prediction/Decisiones/2026-08-10-carlos-elimina-fase3c-confirma-scipion.md`
del vault.

**Codigo eliminado la misma sesion** (todas las secciones de este documento
sobre "Fase A"/clase 1/2/3/Extension 3 arriba quedan como registro
historico de lo que existio y se verifico, no como estado actual):
- `src/structural/fase_a_dispatch.py`, `ddg_estimate.py`,
  `pyrosetta_ptm_patch.py`, `pyrosetta_glycan_patch.py`, `glygen_client.py`,
  `ubiquitin_sumo.py` + sus datos de referencia (`sumo1_reference.pdb`,
  `ubiquitin_reference.pdb`). `src/structural/uniprot_localization_client.py`
  (Fase 3b, via secretora) es la UNICA pieza que sigue viva en ese paquete --
  nunca fue parte de Fase A pese a vivir en el mismo directorio.
- `src/engines/fase_a_engine.py`, `src/engines/_fase_a_runner.py`.
- `src/engines/ptm_annotation.py::select_fase_a_candidates`.
- Bloque `FASE_A_*` completo de `src/config/settings.py` (7 settings:
  `ENABLED`, `TOP_N_PER_TYPE`, `SUPPORTED_PTM_TYPES`, `RESULT_TEMPLATE`,
  `RUNNER_SCRIPT`, `PYTHON_BIN`, `TIMEOUT_SECONDS`).
- `pipeline.py::run_fase_a_pdb_modeling` + el bloque "FASE 3c" del resumen
  en pantalla (`_print_pdb_summary`) + `_FASE_A_RESULT_COLUMNS`.
- Tests dedicados: `test_fase_a_engine.py`, `test_fase_a_dispatch.py`,
  `test_ddg_estimate.py`, `test_pyrosetta_ptm_patch.py`,
  `test_pyrosetta_glycan_patch.py`, `test_glygen_client.py`,
  `test_ubiquitin_sumo.py`. Editados (quitando solo la parte de Fase A, el
  resto del archivo se mantiene): `test_ptm_annotation.py`,
  `test_biological_panel.py`, `test_pipeline_fase1.py`.
- Notebook de Colab (`notebooks/colab_fases_1_3b.ipynb`): quitada la celda
  "Fase 3c (local)" y `FASE_A_ENABLED` de las variables de entorno. El fix
  del mirror de PyRosetta en la celda de DeepPTMPred (2026-08-10, commit
  `23d2dae`) SIGUE siendo necesario -- `DeepPTMPred/predict.py` (Fase 2,
  motor obligatorio) importa `pyrosetta` a nivel de modulo para SASA por
  residuo, sin relacion con Fase A/3c.
- `README.md`: seccion "Fase A" quitada de Arquitectura/Instalacion/
  Licencias, disclaimer de tipo de cadena K48/K63 quitado de "Alcance e
  interpretacion" (dependia de la columna `fase_a_cadena_tipo_aviso`, que ya
  no existe).

**Prioridad inmediata antes de la integracion a Scipion** (pedido explicito
de Enzo, "dejar todo sin friccion"): cachear los 5 entornos conda del
notebook de Colab con `conda-pack` en Drive -- IMPLEMENTADO mas tarde la
misma sesion (commit `58e66dc`, ver seccion "conda-pack" mas abajo) -- y una
auditoria de robustez end-to-end del pipeline actual antes de darlo por
listo para portar a Scipion -- tambien completada la misma sesion
(`/code-review` + pase de seguridad manual + corrida real end-to-end, 6
referencias colgantes a `glygen_client.py` encontradas y corregidas, commit
`8b3473c`).

## StackGlyEmbed ELIMINADO del proyecto (2026-08-10, mismo dia, sesion posterior)

Segunda decision de Carlos en la misma reunion de feedback (ver seccion de
Fase A/3c arriba para la primera): **StackGlyEmbed queda eliminado por
completo del proyecto** -- motivo real: dependia del venv de un proyecto
hermano (`B-Cell-Epitope-Prediction/StackGlyEmbed/.venv-stackglyembed`) +
pesos de un tercero (`scipion-chem-tmbed`, ProtT5 ~3GB), friccion real de
cara a la integracion a Scipion (un protocolo de Scipion no puede depender
de rutas absolutas de OTRO proyecto local). Confirma el mismo patron que ya
motivo la eliminacion de Fase A/3c: friccion de dependencias antes de
portar codigo a Scipion, no un problema de calidad del motor en si.

**Impacto real en el consenso de N-glicosilacion (Camino PDB)**: el
consenso de `n_linked_glycosylation` (decision 2026-08-06) ya estaba
disenado para degradar con cualquier subconjunto de motores disponibles
(`nglyco_consensus_active = pdb_path is not None and Settings.EMNGLY_ENABLED`,
independiente de StackGlyEmbed) -- quitar StackGlyEmbed simplifica el
consenso a 2 motores (DeepMVP+EMNGly) en vez de 3, sin necesidad de tocar
`NGLYCO_CONSENSUS_MIN_ENGINES` (sigue en 2): en la practica ahora exige que
AMBOS motores pasen, en vez de 2 de 3. Verificado leyendo
`_apply_nglyco_consensus` antes de tocar nada -- no fue una suposicion.

**Codigo eliminado**:
- `src/engines/stackglyembed_engine.py`, `src/engines/_stackglyembed_runner.py`.
- `src/engines/ptm_annotation.py`: import de `get_nglyco_corroboration`,
  parametro `enable_stackglyembed` de `annotate_fasta_path`/`annotate_pdb_path`
  (breaking change de API interna, aceptado -- no hay consumidores externos
  del modulo), funcion `_add_stackglyembed_corroboration` completa (rol
  informativo en Camino FASTA, ya no existe en absoluto), rama StackGlyEmbed
  de `_apply_nglyco_consensus` (rol de consenso en Camino PDB).
- Bloque `STACKGLYEMBED_*` completo de `src/config/settings.py` (7 settings:
  `ENABLED`, `PYTHON_BIN`, `RUNNER_SCRIPT`, `MODELS_DIR`, `T5_MODEL_PATH`,
  `ESM_MODEL_NAME`, `TIMEOUT_SECONDS`) -- incluia las 3 rutas absolutas al
  proyecto hermano que causaban la friccion real.
- `pipeline.py`: los 2 call sites que pasaban
  `enable_stackglyembed=Settings.STACKGLYEMBED_ENABLED`.
- Tests: `test_stackglyembed_engine.py` borrado completo; recortado en
  `test_ptm_annotation.py` el bloque entero de corroboracion informativa
  FASTA (8 tests) + 2 tests del "pathway generico" en Camino PDB que ya no
  existe (el fallback a StackGlyEmbed cuando `EMNGLY_ENABLED=False`); los
  tests de consenso de 3 motores reescritos a 2 motores en vez de borrados
  (cobertura equivalente, sin StackGlyEmbed). 338 tests finales (desde 357).
- Notebook de Colab: quitada la Seccion 10 completa ("StackGlyEmbed
  (opcional, desactivado)") y la linea `STACKGLYEMBED_ENABLED` de las
  variables de entorno.
- `README.md`: quitado de Arquitectura (consenso ahora DeepMVP+EMNGly),
  Corroboracion opcional, Alcance e interpretacion, Estado actual, e
  Instalacion (el bloque completo de export de rutas del proyecto hermano).

Commit pendiente (sesion en curso). Notebook probado localmente (sintaxis,
338 tests pasando -- ver arriba), pendiente de confirmar en una corrida real
de Colab.
