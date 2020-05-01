//
// Created by janos on 22.04.20.
//

#include "hash.h"
#include "robin_map.h"

struct A {
    int a,b;
};

struct Hash {
    std::size_t operator()(A& a){
        return hash_int(reinterpret_cast<uint64_t&>(a));
    }
};
int main(){
    tsl::robin_map<A, int, Hash> map;
}
