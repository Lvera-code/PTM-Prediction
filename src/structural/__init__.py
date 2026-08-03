"""Fase A / Extension 3: modelado estructural real de PTMs via PyRosetta.

Los modulos de este paquete NUNCA se importan desde el proceso principal del
pipeline (``pipeline.py``) -- requieren ``pyrosetta``, presente unicamente en
el conda env dedicado ``deepptmpred`` (ver ``STATUS.md``). El unico
consumidor es ``src/engines/_fase_a_runner.py`` (mismo patron que
``_deepptmpred_runner.py``): un script standalone invocado por subprocess
con el interprete de ese conda env, nunca importado directamente por
``src/engines/fase_a_engine.py`` (que corre en el proceso principal, sin
pyrosetta). ``fase_a_dispatch.py`` (dentro de este paquete) es el punto de
entrada que ``_fase_a_runner.py`` importa para rutear cada sitio al modulo
correcto segun su tipo de PTM.
"""
