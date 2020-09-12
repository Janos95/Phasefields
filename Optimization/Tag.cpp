//
// Created by janos on 7/14/20.
//

#include "Tag.h"
#include "Types.h"
#include "StlAlgorithm.h"

#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

static Array<size_t> tags;

size_t getTag() {
    size_t tag = Invalid;
    for(size_t i  = 1; i < tags.size(); ++i, ++tag) {
        if(tags[i] > tags[i - 1] + 1) {
            tag = tags[i - 1] + 1;
            break;
        }
    }

    if(tag == Invalid)
        tag = tags.empty() ? 0 : tags.back() + 1;


    arrayAppend(tags, ++tag);
    std::sort(tags.begin(), tags.end());
    return tag;
}

bool deleteTag(size_t tag) {
    auto it = std::lower_bound(tags.begin(), tags.end(), tag);
    if(it < tags.end() - 1) /* if the tag is not in the last position memmove everything back by one */
        std::memmove(it, it + 1, sizeof(size_t)*(tags.end() - (it + 1)));
    else if(it != tags.end()) /* if the tag is in the last position just pop the back */
        arrayResize(tags, tags.size() - 1);
    else return false;
    return true;
}

}
