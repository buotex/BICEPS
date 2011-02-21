#include "ets.h"






int main(int argc, char *argv[]){
    try{
        Ets(argc,argv);
    }
    catch (std::string & e){
        cout << e << endl;
    }
    catch (const std::exception & e)
    {
        cout << e.what();
    }
    



}

