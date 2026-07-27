# PTM-Prediction

Pipeline de prediccion de zonas de modificacion post-traduccional (PTM) a
partir de FASTA o PDB/mmCIF. Segundo proyecto del CNB (encargado por Carlos
Oscar Sorzano), independiente de
[BCell-Epitope-Prediction](https://github.com/Lvera-code/BCell-Epitope-Prediction).

## Arquitectura

- **Camino FASTA**: Fase 1 (saneamiento) -> **DeepMVP** (motor unico, 6 tipos
  de PTM: fosforilacion, acetilacion, metilacion, sumoilacion,
  ubiquitinacion, N-glicosilacion).
- **Camino PDB**: Fase 1.5 (extraccion de secuencia ATMSEQ + mapeo de
  posiciones via `gemmi`) -> consenso **DeepMVP + DeepPTMPred** (17 tipos de
  PTM; DeepPTMPred exige `pdb_path` obligatorio, sin modo solo-secuencia,
  de ahi la asimetria entre caminos).
- **Fase 3 (nucleo)**: anotacion/filtrado + logica de decision de flujo sobre
  las predicciones de Fase 3a. Diseno cerrado, implementacion pendiente (ver
  STATUS.md).

## Estado actual (2026-07-27)

Fase 1 y Fase 1.5 implementadas y con tests. Fase 3 (motores DeepMVP/
DeepPTMPred + nucleo de anotacion) todavia no construida — ver `STATUS.md`.

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

## Decisiones de arquitectura

Documentadas en el vault, no en este repo:
`01-Proyectos/PTM-Prediction/Decisiones/`.
