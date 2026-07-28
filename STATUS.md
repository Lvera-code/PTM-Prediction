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
- **DeepPTMPred no declara licencia** en su repo (a diferencia de DeepMVP,
  GPL-3.0) — verificar con Carlos antes de cualquier uso mas alla de
  investigacion/TFG. Esto sigue sin resolver.
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
  la pose parcheada relajada localmente (mismo protocolo en ambas). Corrida
  real sobre el mismo sitio: WT relax+score 59s (score 1901.62), MUT
  relax+score 68s (score 1845.40), ddG = **-56.22** (parche estabiliza segun
  esta metrica). **Salvedad real no resuelta**: Tau es intrinsecamente
  desordenada (plDDT medio 49.34 en este modelo AlphaFold, confirmado en el
  log de PyRosetta) y el script hace UN solo relax por estado (sin repetir
  ni promediar) -- el numero es real pero ruidoso/poco fiable como estimado
  de estabilidad en una region de baja confianza estructural; un protocolo
  de produccion querria repetir el relax varias veces por estado y
  promediar/tomar el minimo, no implementado todavia (cada relax ya toma
  ~1 min sobre una proteina de 758 residuos, multiplicar por N repeticiones
  es un coste real a decidir).
- Ninguno de los dos scripts esta conectado a `pipeline.py` todavia --
  sigue la decision 2026-07-27 de que D (el filtro) no rutea a Extension
  3/Fase A porque esas fases no eran parte del nucleo. Se invocan por
  separado, manualmente o desde un futuro orquestador de Fase A/Extension 3.

**Fase A clase 2 (glicosilacion)**: el bloqueante de la documentacion
(403 en `docs.rosettacommons.org`, ver decision del 28-07) sigue sin
resolverse por fetch directo (confirmado de nuevo hoy, mismo 403), PERO se
resolvio por una via mejor: `pyrosetta.rosetta.protocols.carbohydrates.
GlycanTreeModeler` existe y es importable de verdad en el PyRosetta
instalado (introspeccion real, no documentacion) -- confirma que es
utilizable via la API de Python sin coste de instalacion adicional, tal
como se esperaba. Con eso resuelto, `GlycanTreeModeler` es la opcion mas
natural sobre `Glycosylator` (misma logica de "cero coste nuevo" que ya
goberno el resto de Fase A) -- pero la implementacion real (definir el
arbol de glicano a modelar, que ni DeepMVP ni DeepPTMPred predicen -- solo
predicen el SITIO, no la composicion del glicano) no se ha empezado.

**Fase A clase 3 (ubiquitinacion/sumoilacion)**: sin tocar, deliberadamente
la ultima segun el orden decidido, ver razonamiento en la decision del
28-07 (conjugacion de proteina entera, ni Rosetta tiene buena precision
ahi, sin herramienta empaquetada encontrada).

67 tests (`pytest tests/`) siguen pasando tras el parche de
`_deepptmpred_runner.py` -- ningun test nuevo para `src/structural/`
(requiere `pyrosetta`, ausente del entorno donde corre la suite principal,
mismo criterio que el resto de codigo dependiente de motores externos:
verificacion real documentada aqui en vez de mock).

## Proximos pasos reales

1. ~~Decidir el gap de calibracion `fpr` de DeepMVP~~ — RESUELTO 2026-07-28
   noche (ver seccion arriba). `pipeline.py` corre end-to-end de verdad en
   Camino PDB (572 sitios, 98 con consenso real). En una maquina nueva:
   correr `conda run -n deepmvp python scripts/generate_deepmvp_calibration.py`
   una vez despues de instalar los pesos de DeepMVP, antes del primer uso
   real del pipeline.
2. Repetir/promediar el relax en `ddg_estimate.py` (varias corridas por
   estado) si se decide usarlo como metrica real de produccion, no solo de
   verificacion puntual.
3. Fase A clase 2 (glicosilacion): definir de donde sale la composicion del
   glicano a modelar con `GlycanTreeModeler` (ninguno de los dos motores la
   predice) antes de poder implementar de verdad.
4. Fase A clase 3 (ubiquitinacion/sumoilacion): sin empezar, deliberadamente
   al final.
5. Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion —
   deliberadamente sin integrar por ahora (decision 2026-07-26).
