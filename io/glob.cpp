//
// Created by jmeny on 08.11.19.
//

#include "glob.hpp"

#include <Corrade/Utility/Directory.h>

#include <string_view>
#include <string>

std::vector<std::string> glob(
        const std::string& path,
        const std::string& extension,
        std::vector<std::string> paths,
        bool recurse)
{
    using namespace Corrade::Utility;
    paths.clear();
    for(const auto& path: Directory::list(path))
    {
        if(Directory::isDirectory(path) && recurse){
            glob(path, extension, paths, recurse);
            continue;
        }
        auto [_, ext] = Directory::splitExtension(path);

        if(extension == ext)
            paths.push_back(std::move(path));
    }

    return paths;
}