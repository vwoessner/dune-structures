add_executable(structured2d_p2 EXCLUDE_FROM_ALL structured2d_p2.cc)
add_executable(structured3d_p1 EXCLUDE_FROM_ALL structured3d_p1.cc)

add_executable(cell3d_p1 EXCLUDE_FROM_ALL cell3d_p1.cc)
add_executable(cell3d_p2 EXCLUDE_FROM_ALL cell3d_p2.cc)

# Register all these apps for the frontend target
blocklab_add_app_to_frontend(TARGET structured2d_p2)
blocklab_add_app_to_frontend(TARGET structured3d_p1)
blocklab_add_app_to_frontend(TARGET cell3d_p1)
blocklab_add_app_to_frontend(TARGET cell3d_p2)
