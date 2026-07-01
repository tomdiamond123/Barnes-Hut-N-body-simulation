Barnes-Hut N body simulation using C++ and SFML.
Simulates 10,000 bodies at ~60 fps.

Required - Install CMake:
https://cmake.org/download/

Make sure to compile the release build for full performance benifits
- CMake: Select Variant command in VSCode
- Or make sure -O3 flag is present when building from command line with gcc

If using GCC/G++:

Required - install tbb for parallel execution:

Run: 

pacman -S mingw-w64-ucrt-x86_64-tbb

in the MSYS2 UCRT64 Terminal

If this doesn't work or you get errors, change this code in the main game loop, this will result in lower frame rates:

std::for_each(std::execution::par_unseq, bodies.begin(), bodies.end(), 
                        [quadtree](Body& body){
                            if (!body.active) return;

                            quadtree->updateBodyPosition(body);

                            if (body.position.first < 0 || body.position.second < 0 || 
                                body.position.first > SIMSIZE || body.position.second > SIMSIZE){
                                    body.active = false;
                                    return;
                                }
                        });

        for (Body& body: bodies){
            if (!body.active) continue;
            bodyShape.setPosition({static_cast<float>(body.position.first * SCALE), static_cast<float>(body.position.second * SCALE)});
            window.draw(bodyShape);
        }

to this:

for (Body& body: bodies){
            if (!body.active) continue;

            if (body.position.first < 0 || body.position.second < 0 || 
                body.position.first > SIMSIZE || body.position.second > SIMSIZE){
                    body.active = false;
                    continue;
                }

            quadtree->updateBodyPosition(body);
            bodyShape.setPosition({static_cast<float>(body.position.first * SCALE), static_cast<float>(body.position.second * SCALE)});
            window.draw(bodyShape);
        }

And remove tbb12 from CMakeLists.txt on this line:

target_link_libraries(barneshut PRIVATE SFML::Graphics tbb12)


Barnes-Hut Simulation Explanation:
https://arborjs.org/docs/barnes-hut

SFML Documentation:
https://www.sfml-dev.org/documentation/3.1.0/
https://www.sfml-dev.org/tutorials/3.1/getting-started/cmake/

SFML Template:
https://github.com/SFML/cmake-sfml-project

Set Up CMake on VSCode - tutorial for linux but works on all platforms
https://code.visualstudio.com/docs/cpp/cmake-linux

Timing Class:
https://www.learncpp.com/cpp-tutorial/timing-your-code/