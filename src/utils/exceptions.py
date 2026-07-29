"""Jerarquia de excepciones especificas del pipeline para control de flujo."""


class PipelineError(Exception):
    """Clase base para todos los errores controlados del pipeline."""


class InvalidSequenceError(PipelineError):
    """Una secuencia FASTA individual tiene residuos invalidos o longitud insuficiente.

    Recuperable a nivel de registro: el saneamiento lo captura internamente y
    descarta unicamente la secuencia afectada, sin detener el lote. Fatal a
    nivel de archivo si NINGUN registro sobrevive el saneamiento.
    """


class FastaFormatError(PipelineError):
    """El archivo FASTA de entrada no cumple la sintaxis minima (sin cabeceras '>')
    o contiene accessions duplicados. Fatal: detiene el pipeline antes de
    iniciar cualquier fase.
    """


class StructureParsingError(PipelineError):
    """Fallo al parsear una estructura de entrada (PDB/mmCIF) en Fase 1.5.

    Cubre archivos corruptos o sin sintaxis valida, ausencia de cualquier
    cadena proteica con al menos un residuo de aminoacido valido en la
    estructura, y fallos al resolver la estrategia de seleccion de cadena
    configurada (``Settings.PDB_CHAIN_SELECTION_STRATEGY``). Es fatal: no hay
    manera segura de continuar sin una cadena de referencia valida.
    """


class InputRoutingError(PipelineError):
    """Fallo al determinar el tipo de un archivo de entrada (FASTA vs estructura).

    Cubre un archivo cuyo tipo no puede determinarse con confianza ni por
    extension ni por contenido. Es fatal: detiene el pipeline antes de correr
    cualquier fase.
    """


class EngineExecutionError(PipelineError):
    """Error durante la ejecucion de un motor de prediccion de PTMs (subprocess
    o forward pass). Recuperable a nivel de lote: se loggea y se propaga para
    que el orquestador decida si continuar con los lotes restantes o abortar.
    """


class DeepMVPExecutionError(EngineExecutionError):
    """Fallo al ejecutar DeepMVP localmente (motor unico Camino FASTA, motor 1
    de 2 en consenso Camino PDB).
    """


class DeepPTMPredExecutionError(EngineExecutionError):
    """Fallo al ejecutar DeepPTMPred localmente (motor 2 de 2 en consenso,
    Camino PDB unicamente: requiere pdb_path obligatorio, sin modo
    solo-secuencia).
    """
