# STATUS

Estado actual del proyecto, no un diario de sesiones — reescribir en vez de
acumular. Ver `01-Proyectos/PTM-Prediction/` en el vault para el historial
de decisiones.

## Hecho

- **Fase 1** (`src/utils/fasta_parser.py`): saneamiento de FASTA crudo
  (mayusculas, rechazo fatal de residuos no canonicos, deteccion de
  accessions duplicados). Reutiliza el patron de
  `BCell-Epitope-Prediction/src/utils/fasta_parser.py`, con una diferencia
  real: la politica de rechazo de residuos no canonicos aqui es un default
  CONSERVADOR, no una verificacion confirmada del comportamiento real de
  DeepMVP (a diferencia de BepiPred-3.0 en proyecto 1, confirmado por
  lectura directa del codigo). Revisar cuando se construya
  `deepmvp_engine.py`.
- **Fase 1.5** (`src/utils/structure_parser.py`): extraccion de secuencia
  ATMSEQ + mapeo de posiciones desde PDB/mmCIF via `gemmi`. Identico al de
  proyecto 1 (logica 100% generica, no especifica de epitopos).
- **Enrutador de input** (`src/utils/input_router.py`): identico al de
  proyecto 1 (extension + sniffing de contenido).
- **Orquestador** (`pipeline.py`): CLI minima que corre Fase 1 o 1.5 segun
  el tipo de input detectado. Fase 3 todavia no esta enganchada (imprime un
  aviso, no falla en silencio).
- 34 tests (`pytest tests/`, sin binarios/modelos externos).
- Repo local (`git init`), sin remoto todavia — decidido 2026-07-27, se
  crea en GitHub (`Lvera-code/PTM-Prediction`, publico) cuando haya algo
  funcional end-to-end.

## Diseno cerrado, pendiente de construir

**Fase 3 (nucleo)** — diseno completo en
`01-Proyectos/PTM-Prediction/Decisiones/2026-07-27-diseno-nucleo-fase3-anotacion-flujo.md`:

- Esquema de salida: `accesion | posicion | residuo_wt | tipo_PTM | motor(es)
  | score_DeepMVP | score_DeepPTMPred (nullable) | consenso (bool, solo
  Camino PDB) | ventana (nullable) | camino (FASTA/PDB)`.
- Umbral por herramienta (no global), score crudo siempre conservado.
- Los ~11 tipos exclusivos de DeepPTMPred en Camino PDB se incluyen en el
  nucleo, marcados `consenso=false`.
- Logica de flujo (D): filtro/prioridad de responsabilidad unica (umbral
  generico pasa/no-pasa), sin rutas a Extension 3 (ΔΔG) ni Fase A
  (modelado estructural) — esas fases no existen todavia.

**Siguiente paso real**: construir `src/engines/deepmvp_engine.py` — primero
verificar por lectura directa del repo (`raw.githubusercontent.com`, no
resumen de buscador) el CLI real, el umbral recomendado publicado, y la
tolerancia real a residuos no canonicos (para ajustar o confirmar la
politica conservadora de `fasta_parser.py` de hoy). Repo/pesos de DeepMVP y
DeepPTMPred ya verificados como instalables el 2026-07-26 (ver decision de
esa fecha), pero ninguno de los dos se ha clonado/instalado todavia en esta
maquina.

## Explicitamente fuera del nucleo (no descartado, agendado)

- Extension 3 (ΔΔG / impacto de estabilidad).
- Fase A (modelado estructural real del PTM, Camino PDB unicamente).
- Cross-validacion con StackGlyEmbed (proyecto 1) para N-glicosilacion —
  ambos motores quedan independientes por ahora.
