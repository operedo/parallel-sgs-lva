#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" {
#endif
//    class Foo;
//    typedef Foo FOO;
//#else
//    // From the C side, we use an opaque pointer.
//    typedef struct FOO FOO;
//#endif
//
//// Constructor
//FOO* create_foo(int a, int b);
//
//// Destructor
//void delete_foo(FOO* foo);
//
//// The const qualificators maps from the member function to pointers to the
//// class instances.
//int foo_bar(const FOO* foo, int c);
//double foo_baz(const FOO* foo, double d);
//
//void foo_speaker(const char* s);

//void dijkstra(const char* s);
void dijkstra(int* n, int* a, int* b, int* c, double* d, int* e, int* f, double* g);

#ifdef __cplusplus
}
#endif
