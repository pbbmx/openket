import sys
import os

# Añade el directorio que contiene 'openket' al sys.path
sys.path.append(os.path.abspath('/home/ultrxvioletx/GitHub'))
print(sys.path)

# Ahora intenta la importación
try:
    from openket.core.diracobject import DiracObject, Ket, Bra, Operator, CreationOperator, AnnihilationOperator
    print("¡Importación Exitosa!")

    # Crea instancias para verificar que las clases están disponibles
    ket = Ket()
    print(ket)

except ImportError as e:
    print(f"Error al importar: {e}")