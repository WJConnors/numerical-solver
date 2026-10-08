add_executable(laplace
	main.cpp
)

target_link_libraries(laplace PRIVATE
	numerical-common
)
