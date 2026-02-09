# L'outil CMake permet de générér automatiquement le fichier Makefile d'un projet. Pour utiliser cet outil, il faut écrire un fichier nommé CMakeLists.txt. L'intérêt est que son écriture est relativement simple.
# Où copier le fichier ?
#        Il faut le placer à la racine du projet (à l'endroit où se trouve les dossiers src et data)
#
#    Comment compiler la première fois ?
#       1. A la racine de votre projet (où vous avez placé votre fichier CMakeLists.txt), créer un dossier appelé build avec la commande suivante : "mkdir build")
#       2. Se placer dans le dossier build avec la commande suivante : "cd build"
#       3. Créer le fichier Makefile associé à votre projet avec la commande suivante : "cmake .." (compilation en mode RELEASE)  OU "cmake -DCMAKE_BUILD_TYPE=Debug .." (compilation en mode DEBUG) OU "cmake -DCMAKE_BUILD_TYPE=Release .." (compilation en mode RELEASE)
#       4. Si la création du fichier Makefile a fonctionné (vous voyez dans la console "-- Generation done"), vous pouvez simplement lancer la compilation comme vous le faites habituellement. Pusique la cible s'appelle ici gurobiSolverCWLP (il est possible de modifier ce nom dans le fichier CMakeLists.txt), il faut taper la commande suivante : make gurobiSolverCWLP.
#       5. Si la compilation n'échoue pas, un exécutable appelé gurobiSolverCWLP a été créé. Vous pouvez le lancer en tapant la commande suivante :  ./gurobiSolverCWLP

#    Comment compiler après avoir effectué des modifications ?
#        si vous n'avez pas ajouté pas de nouveaux fichiers au dossier src,il suffit de taper à nouveau la commande suivante : make gurobiSolverCWLP (en se placant d'abord dans le dossier build)
 #       si vous avez ajouté de nouveaux fichiers au dossier src,il suffit de reprendre à partir de l'étape 2 ci-dessus

cmake_minimum_required(VERSION 3.14)
project(TSP_Gurobi)

set(CMAKE_CXX_STANDARD 14)

#------------------- CREATION DES VARIABLES POUR LA BIBLIOTHEQUE GUROBI -----------------------
message(STATUS "CMAKE_SYSTEM_NAME : ${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux") #Pour Linux
	set(GUROBI_ROOT "/opt/gurobi1301/linux64")  #A Modifier si votre chemin est différent
endif()	
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows") #Pour Windows
	set(GUROBI_ROOT "C:/gurobi952/win64") #A Modifier si votre chemin est différent
endif()
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin") #Pour MacOS
	set(GUROBI_ROOT "/Library/gurobi952/mac64") #A Modifier si votre chemin est différent
endif()
	
message(STATUS "GUROBI_ROOT : ${GUROBI_ROOT}")

# On créé une variable stockant le dossier où se situent les fichiers en-tête du solveur GUROBI (pour la compilation)
set(GUROBI_INCLUDE_DIR "${GUROBI_ROOT}/include")
message(STATUS "GUROBI_INCLUDE_DIR : ${GUROBI_INCLUDE_DIR}") #On affiche le dossier précédent (pour DEBUG)

# On créé des variables stockant l'emplacement des biblothèques/librairies du solveur GUROBI (pour l'édition des liens)
# ATTENTION -> Selon la version utilisée, il faut modifier le nom de la librairie ci-dessous
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux") #Pour Linux
	set(GUROBI_LIBRARY_1 "${GUROBI_ROOT}/lib/libgurobi_c++.a")
	set(GUROBI_LIBRARY_2 "${GUROBI_ROOT}/lib/libgurobi130.so")
endif()
if (${CMAKE_SYSTEM_NAME} MATCHES "Windows") #Pour Windows
	set(GUROBI_LIBRARY_1 "${GUROBI_ROOT}/lib/gurobi95.lib")
	if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    	set(GUROBI_LIBRARY_2 "${GUROBI_ROOT}/lib/gurobi_c++mdd2019.lib")
	else()
    	set(GUROBI_LIBRARY_2 "${GUROBI_ROOT}/lib/gurobi_c++md2019.lib")
	endif()
endif()
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin") #Pour MacOS
	set(GUROBI_LIBRARY_1 "${GUROBI_ROOT}/lib/libgurobi_c++.a")
	set(GUROBI_LIBRARY_2 "${GUROBI_ROOT}/lib/libgurobi91.dylib")
endif()
set(GUROBI_LIBRARIES ${GUROBI_LIBRARY_1} ${GUROBI_LIBRARY_2})
message(STATUS "GUROBI_LIBRARIES : ${GUROBI_LIBRARIES}") #On affiche les chemins vers les biblothèques (pour DEBUG)

#------------------- PARTIE PRINCIPALE -----------------------
# On indique où trouver les définitions des fonctions du solveur GUROBI
include_directories(${GUROBI_INCLUDE_DIR})
include_directories(include)

# Recherche tous les fichiers .cpp dans le dossier src (stocke tout cela dans la variable SRCS)
file(GLOB_RECURSE SRCS src/*.cpp)
message(STATUS "SRCS : ${SRCS}") #Affiche tous les fichiers cpp trouvés (pour DEBUG)

# On indique que l'on veut un exécutable "gurobiSolverCWLP" compilé à partir des fichiers stockés dans SRCS
add_executable(TSP_Gurobi ${SRCS})

# On indique où trouver les biblothèques du solveur GUROBI
target_link_libraries(TSP_Gurobi ${GUROBI_LIBRARIES})

# On ajoute un lien symbolique vers le dossier data dans le dossier où se situera l'exécutable
# (cela permet d'indiquer le chemin d'une instance contenu dans le dossier data en indiquant simplement un chemin relatif et non absolu)
if (UNIX)  #Pour MacOS et Linux
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data/ ${CMAKE_CURRENT_BINARY_DIR}/data)
else()  #Pour Windows
    file(TO_NATIVE_PATH "${CMAKE_CURRENT_BINARY_DIR}/data" _dstDir)
    file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/" _srcDir)
    execute_process(COMMAND cmd.exe /c mklink /J "${_dstDir}" "${_srcDir}")
endif()
