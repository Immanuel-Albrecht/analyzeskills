

FILE(GLOB SrcSubdirs ./boolean/*.c )


add_executable(../analyzeskills analyzeskills.c ${SrcSubdirs})

target_link_libraries(../analyzeskills
            cprops )

