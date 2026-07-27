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

## Riesgos de instalacion identificados, no verificados todavia

- **DeepPTMPred usa `tensorflow-addons`** (para su loss function
  `SigmoidFocalCrossEntropy`), paquete archivado/deprecado por Google desde
  2024, contra TensorFlow 2.15 — no se ha intentado instalar todavia, riesgo
  real de incompatibilidad no descartado.
- **DeepPTMPred requiere PyRosetta** (licencia academica gratuita via
  registro, mismo matiz que otras herramientas de licencia restringida en
  proyecto 1), checkpoint ESM-2 de ~2.5GB, GPU recomendada, ~50GB de disco.
  Instalacion no iniciada en esta maquina.
- **DeepPTMPred no declara licencia** en su repo (a diferencia de DeepMVP,
  GPL-3.0) — verificar con Carlos antes de cualquier uso mas alla de
  investigacion/TFG.
- El runner propio (`_deepptmpred_runner.py`) esta escrito y con la logica
  verificada contra el codigo fuente real, pero NUNCA se ha ejecutado
  contra el entorno real (sin PyRosetta/TF/fair-esm instalados aqui) — solo
  probado con `subprocess.run` mockeado.

## Proximos pasos reales

1. Clonar DeepMVP + descargar sus pesos, clonar DeepPTMPred + descargar
   checkpoint ESM-2 + instalar PyRosetta, verificar ambos venvs dedicados
   (stacks incompatibles entre si: DeepMVP Python 3.7/TF 2.4, DeepPTMPred
   Python 3.10/TF 2.15).
2. Correr el pipeline real end-to-end sobre un caso real (validar que el
   runner de DeepPTMPred funciona de verdad, no solo mockeado).
3. Extension 3 (ΔΔG) y Fase A (modelado estructural real) — diseno cerrado
   el 26-07, implementacion no empezada, deliberadamente pospuestas.
4. Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion —
   deliberadamente sin integrar por ahora (decision 2026-07-26).
