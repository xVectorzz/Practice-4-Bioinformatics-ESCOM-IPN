# Práctica 4: Simulación de Fuerzas 2D (Dinámica Molecular)

Esta práctica realiza una simulación física de interacciones entre moléculas de agua y un ion de sodio usando fuerzas de Hooke, Lennard-Jones y Coulomb, visualizado en tiempo real con OpenCV.

## Instrucciones de uso

Asegúrate de tener `Practica4.py` en tu carpeta. Abre tu terminal y corre estos comandos en orden para preparar el entorno, instalar lo necesario y ejecutar el script:

**En Ubuntu (WSL / Linux):**
```bash
# Crear entorno, activar, instalar dependencias y ejecutar
python3 -m venv bioenv
source bioenv/bin/activate
pip install numpy opencv-python periodictable
python3 Practica4.py
```

**En Windows (CMD):**
```cmd
# Crear entorno, activar, instalar dependencias y ejecutar
python -m venv bioenv
bioenv\Scripts\activate
pip install numpy opencv-python periodictable
python Practica4.py
```
