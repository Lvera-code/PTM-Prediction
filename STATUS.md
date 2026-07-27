# STATUS

Estado actual del proyecto, no un diario de sesiones — reescribir en vez de
acumular. Ver `01-Proyectos/PTM-Prediction/` en el vault para el historial
de decisiones.

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

- **DeepMVP — INSTALADO Y FUNCIONAL** (sin pesos, esos si requieren
  descarga manual del Shiny app). Entorno conda real `deepmvp`
  (`/home/enzo/miniconda3/envs/deepmvp/bin/python`, Python 3.7.10):
  `pip install -r requirements.txt` completo sin errores (TensorFlow
  2.4.2 CPU, numpy 1.19.5, todo el stack). `python DeepMVP.py predict -h`
  imprime exactamente el CLI documentado arriba. Verificado tambien el
  caso de error real: con `-m` apuntando a una carpeta de modelos VACIA,
  DeepMVP.py revienta con `ValueError: No objects to concatenate` (traza
  cruda de pandas, opaca) — confirma que el chequeo proactivo de
  `DeepMVPEngine._validate_installation` (falla antes, con mensaje
  accionable) es necesario, no cosmetico. Un warning benigno de numpy/API
  mismatch aparece en stderr pero no afecta la ejecucion.
- **DeepPTMPred — dependencias mas riesgosas YA VERIFICADAS, instalacion
  completa NO intentada** (PyRosetta requiere su propio instalador +
  registro, ESM-2 checkpoint ~2.5GB de descarga manual, ~50GB disco — fuera
  del tiempo disponible hoy). Verificado en un venv Python 3.10 aislado:
  `tensorflow==2.15` + `tensorflow-addons==0.23.0` instalan Y IMPORTAN sin
  error (`from tensorflow_addons.losses import SigmoidFocalCrossEntropy`
  funciona) — solo un `UserWarning` de deprecacion esperado (TFA en modo
  mantenimiento desde 2024, fin de vida anunciado mayo 2024, pero el
  paquete sigue siendo instalable e importable). Esto DESCARTA el riesgo
  de incompatibilidad que se habia flageado sin verificar. `fair-esm` +
  `torch` tambien instalan e importan limpio. **Unico bloqueante real
  restante: PyRosetta.** Se intento la instalacion automatica real via
  `pip install pyrosetta-installer` + `install_pyrosetta()` (el metodo que
  documenta el propio README de DeepPTMPred) — el paquete instalador SI se
  instala, pero la descarga del wheel real falla en esta maquina: el mirror
  por defecto (`west.rosettacommons.org`) responde `404` en la ruta que el
  instalador consulta (`.../latest.html`), y el mirror alternativo
  (`graylab.jhu.edu`) falla la verificacion TLS (cadena de certificados no
  confiable en este entorno). No investigado mas a fondo (puede ser un
  cambio de ruta del lado de PyRosetta, o una restriccion de red/CA de este
  entorno sandboxeado) — instalacion manual necesaria, confirmado que NO es
  automatizable tal cual en esta maquina, no es simplemente "no intentado".
- **DeepPTMPred no declara licencia** en su repo (a diferencia de DeepMVP,
  GPL-3.0) — verificar con Carlos antes de cualquier uso mas alla de
  investigacion/TFG. Esto sigue sin resolver.
- **Confirmado 100% local en toda la cadena (2026-07-27)**: verificado
  leyendo `esm/pretrained.py` de github.com/facebookresearch/esm
  directamente. `pretrained.load_model_and_alphabet(path)` solo llama a
  red (`dl.fbaipublicfiles.com`) si el argumento NO termina en `.pt` (rama
  hub); como el runner siempre pasa una ruta `.pt` local, siempre entra por
  la rama `load_model_and_alphabet_local` (`torch.load()` puro sobre
  disco). Detalle real (no de red, de archivos): esa rama tambien exige un
  companero `<checkpoint>-contact-regression.pt` en el mismo directorio
  (heuristica de fair-esm que no excluye `esm2_*`) — si falta, falla con
  `FileNotFoundError` local al descargar solo el checkpoint principal.
  Verificar al bajar el checkpoint que ese archivo companero tambien este
  disponible junto a el.
- El runner propio (`_deepptmpred_runner.py`) sigue sin ejecutarse contra
  el entorno real completo (falta PyRosetta + checkpoint ESM-2 + pesos) —
  solo probado con `subprocess.run` mockeado en tests.

## Proximos pasos reales

1. Descargar pesos de DeepMVP (manual, `https://deepmvp.ptmax.org/`,
   Shiny app) y apuntar `DEEPMVP_MODEL_DIR` — el resto de la instalacion
   ya esta lista y verificada (env conda `deepmvp` funcional).
2. Instalar PyRosetta manualmente (la instalacion automatica via
   `pyrosetta-installer` falla en esta maquina, ver arriba — probar desde
   una red/maquina distinta, o instalar via conda con canal academico
   registrado) + descargar checkpoint ESM-2 (~2.5GB) para DeepPTMPred —
   unico bloqueante real restante de ese lado, ya que TF/tensorflow-addons/
   fair-esm/torch se verificaron limpios.
3. Correr el pipeline real end-to-end sobre un caso real una vez ambos
   esten completos (validar que el runner de DeepPTMPred funciona de
   verdad, no solo mockeado).
4. Extension 3 (ΔΔG) y Fase A (modelado estructural real) — diseno cerrado
   el 26-07, implementacion no empezada, deliberadamente pospuestas.
5. Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion —
   deliberadamente sin integrar por ahora (decision 2026-07-26).
