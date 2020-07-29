//
// Created by janos on 7/14/20.
//

#include "tag.h"

#include <set>

static std::set<int> usedTags;

int getTag() {
    int tag = -1;
    auto it = usedTags.begin();
    while(it != usedTags.end()){
        if(*it > tag + 1){
            usedTags.insert(++tag);
            return tag;
        }
        tag = *it;
    }
    usedTags.insert(++tag);
    return tag;
}

void deleteTag(int tag) {
    usedTags.erase(tag);
}
