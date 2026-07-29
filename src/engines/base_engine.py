"""Interfaz abstracta comun para los motores de prediccion de PTMs (DeepMVP, DeepPTMPred)."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, Sequence, TypeVar

TIn = TypeVar("TIn")
TOut = TypeVar("TOut")


class BaseEngine(ABC, Generic[TIn, TOut]):
    """Contrato minimo que debe cumplir cualquier motor de prediccion del pipeline.

    Cada motor (DeepMVP, DeepPTMPred) implementa esta interfaz para permitir
    composicion y sustitucion transparente (patron Strategy), sin acoplar el
    orquestador (``pipeline.py``) a una implementacion concreta.
    """

    @abstractmethod
    def run(self, items: Sequence[TIn], output_dir: Optional[Path] = None) -> List[TOut]:
        """Ejecuta inferencia sobre un lote de datos de entrada.

        Args:
            items: Secuencia homogenea de elementos de entrada para este motor.
            output_dir: Carpeta base donde el motor persiste sus artefactos
                crudos (subcarpeta por accession/stem). Si es ``None``, cada
                motor concreto usa su propio default (``Settings.FASTA_OUTPUT_DIR``
                en DeepMVPEngine/DeepPTMPredEngine, ver src/config/settings.py).
                Parte del contrato desde el principio (audit 2026-07-28: antes
                era un parametro no declarado que ambos motores concretos ya
                requerian en la practica, rompiendo el polimorfismo declarado).

        Returns:
            Lista de resultados, en el mismo orden logico que ``items`` salvo
            que el motor documente explicitamente lo contrario.
        """
        raise NotImplementedError
