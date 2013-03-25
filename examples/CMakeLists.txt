include_directories(${Moirai_INCLUDE_DIRECTORIES})


add_executable (test_pde_mol test_pde_mol.cpp) 
target_link_libraries (test_pde_mol Moirai)

