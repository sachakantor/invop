Linkeados con CPLEX
===================
CVRP y Food_Manufacture han sido implementados en C++, bajos sendos proyectos de CMake, y linkeados contra la librería de CPLEX.

Para compilar dichos códigos:
1) Configurar la ruta a la suite CPLEx local en el archivo CMakeLists.txt . Buscar la variable CPLEX_STUDIO_DIR y definirla para los datos locales.
2) Desde una consola, posicionarse en el directorio raíz del repo
3) Ejecutar el siguiente comando:

mkdir build && cd build && cmake -G "Unix Makefiles" .. && make all

4) Si todo ejecutó bien, los binarios ejecutables correspondientes se encontrarán en

CVRP/bin/cvrp
Food_Manufacture/bin/food_manufacture

5) El formato de ejecución de ambos binarios es

CVRP/bin/cvrp ARCHIVO_INPUT1 [ARCHIVO_INPUT2...]

El ARCHIVO_INPUT es la instancia del problema, que tambien se provee como parte de este archivo. Los dos casos de ejecución son:

cvrp/bin/cvrp cvrp/problem_instance.in
Food_Manufacture/bin/food_manufacture food_manufacture/problem_instance_input_distances.in

