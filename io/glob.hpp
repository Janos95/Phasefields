//
// Created by jmeny on 08.11.19.
//

#pragma once

#include <vector>
#include <string>


std::vector<std::string> glob(
        const std::string& path,
        const std::string& extension,
        std::vector<std::string> paths,
        bool recurse = false);
