#include "Biceps.h"






int main(int argc, char *argv[]){
    try{
        Biceps(argc,argv);
    }
    catch (std::string & e){
        cout << e << endl;
    }
    catch (const std::exception & e)
    {
        cout << e.what();
    }
    



}

