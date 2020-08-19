//
// Created by janos on 7/14/20.
//

#include "Tag.h"

#include <Corrade/Containers/GrowableArray.h>
#include <algorithm>

using namespace Corrade;

static Containers::Array<int> g_tags;

int getTag() {
    int tag = -1;
    std::size_t i = 0;
    for(; i < g_tags.size(); ++i, ++tag) {
        if(g_tags[i] > tag + 1){
            ++tag;
            break;
        }
    }

    // @TODO this can be done more efficiently
    Containers::arrayAppend(g_tags, ++tag);
    std::sort(g_tags.begin(), g_tags.end());
    return tag;
}

void deleteTag(int tag) {
    auto it = std::lower_bound(g_tags.begin(), g_tags.end(), tag);
    if(it < g_tags.end() - 1){
        std::memmove(it, it + 1, sizeof(int)*(g_tags.end() - (it + 1)));
    }
    if(it != g_tags.end()){
        Containers::arrayResize(g_tags, g_tags.size() - 1);
    }
}
