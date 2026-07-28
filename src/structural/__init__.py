"""Fase A / Extension 3: modelado estructural real de PTMs via PyRosetta.

Los modulos de este paquete NUNCA se importan desde el paquete ``src``
principal (mismo patron que ``src/engines/_deepptmpred_runner.py``):
requieren ``pyrosetta``, presente unicamente en el conda env dedicado
``deepptmpred`` (ver ``STATUS.md``). Se invocan siempre por subprocess o
directamente con el interprete de ese env.
"""
