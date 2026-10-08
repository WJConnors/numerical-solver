add_executable(poisson
	main.cpp
)

target_link_libraries(poisson PRIVATE
	numerical-common
)
