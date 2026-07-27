"""Config raiz de pytest: garantiza que la raiz del repo este en sys.path
para que los tests puedan importar 'pipeline' y 'src.*' sin instalar el
paquete. No requiere ningun binario/modelo externo (DeepMVP/DeepPTMPred):
los tests en tests/ cubren unicamente la logica pura de Fase 1/1.5 (parseo,
saneamiento, extraccion de estructura), nunca invocan un motor real.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
