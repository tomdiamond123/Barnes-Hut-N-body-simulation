#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <execution>

using coords = std::pair<double, double>;

// Explanation of Barnes-Hut algorithm - https://arborjs.org/docs/barnes-hut

// Code tagged //toggle is for logging and measuring performance of different parts of the program so can be turned on or off

// Code in --------- Change these --------- blocks can be changed to influence the behaviour of the sim.

// --------- Change these ---------
const int WINDOWSIZE = 1500;
const double STARTINGSIZE {1*(3e8*365.25*24*60*60)}; // size in metres

// buffer round the edge of STARTING SIZE to allow for movement of particles outside of starting area
const double SIMSIZE {STARTINGSIZE*1.5}; 

const int NUMOFBODIES {5000};

const double theta = 0.5; // higher number = less acurate prediction but higher performance, 0.0 = brute force

const double TIMESTEP = 3600*24*5000; // amount of time that passes each frame in seconds

const double STARTMASS = 1e30; // in kilograms

const double EPSILON = SIMSIZE * 0.001; // softens forces produced by objects being very close to eachother
// const double EPSILON = 0;

const sf::Color COLOUR (0, 255, 255);
// --------------------------------

const double SCALE {WINDOWSIZE/SIMSIZE};
const double G {6.67438e-11};
const double PI {3.14159265358979};
// const double STARTINGVELOCITYRANGE = 1e5;

// timing class - https://www.learncpp.com/cpp-tutorial/timing-your-code/
class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<double, std::ratio<1> >;

	std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
	void reset()
	{
		m_beg = Clock::now();
	}

	double elapsed() const
	{
		return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
	}
};

struct Body{
    int id{};
    coords position{};
    double mass{};
    double xVel{};
    double yVel{};
    bool active{true};
};


class Node {
    public:
        coords centreOfMass{};
        double mass{};
        double length{};

        //coords of top left and bottom right of node
        coords boundTopLeft{};
        coords boundBottomRight{};

        bool isLeaf{};

        Node* topLeft;
        Node* topRight;
        Node* bottomLeft;
        Node* bottomRight;

        int bodyId{-1};

        // constructor
        Node(int id, coords centreOfMass, double mass, double length, coords boundTopLeft, coords boundBottomRight, bool isLeaf) 
        : bodyId(id), centreOfMass(centreOfMass), mass(mass), length(length), boundTopLeft(boundTopLeft), boundBottomRight(boundBottomRight), 
        isLeaf(isLeaf), topLeft(nullptr), topRight(nullptr), bottomLeft(nullptr), bottomRight(nullptr) {}

        //destructor
        ~Node(){
            delete topLeft;
            delete topRight;
            delete bottomLeft;
            delete bottomRight;
        }
};


class QuadTree {
private:
    //pointer to the root of the tree
    Node* root {nullptr};
    const std::vector<Body>& bodies;

    void updateMassAndCentreofMass(const Body& body, Node* node){
        auto [nodeX, nodeY] = node->centreOfMass;
        double mass = node->mass + body.mass;
        nodeX = (nodeX*node->mass + body.position.first*body.mass)/mass;
        nodeY = (nodeY*node->mass + body.position.second*body.mass)/mass;
        node->mass += body.mass;
        node->centreOfMass = {nodeX, nodeY};
    }

    void subdivide(const Body& body, Node* node){
        node->isLeaf = false;
        int oldId = node->bodyId;
        node->bodyId = -1;
        recursiveInsertBody(bodies[oldId], node);
        recursiveInsertBody(body, node);
        updateMassAndCentreofMass(body, node);
    }

    void recursiveInsertBody(const Body& body, Node* root){
        auto [x, y] = body.position;

        double midX = (root->boundTopLeft.first + root->boundBottomRight.first) / 2.0;
        double midY = (root->boundTopLeft.second + root->boundBottomRight.second) / 2.0;

        //Check top left
        if (x <= midX && y <= midY){
            if (root->topLeft == nullptr){
                coords boundBottomRight = {root->boundBottomRight.first-root->length/2, root->boundBottomRight.second-root->length/2};
                root->topLeft = new Node(body.id, body.position, body.mass, root->length/2, root->boundTopLeft, boundBottomRight, true);
            }
            else{
                Node* node = root->topLeft;
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
                    updateMassAndCentreofMass(body, node);
                    recursiveInsertBody(body, node);
                }
            }
        }
        //Check top right
        else if (x > midX && y <= midY){
            if (root->topRight == nullptr){
                coords boundTopLeft = {root->boundTopLeft.first+root->length/2, root->boundTopLeft.second};
                coords boundBottomRight = {root->boundBottomRight.first, root->boundBottomRight.second-root->length/2};
                root->topRight = new Node(body.id, body.position, body.mass, root->length/2, boundTopLeft, boundBottomRight, true);
            }
            else{
                Node* node = root->topRight;
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
                    updateMassAndCentreofMass(body, node);
                    recursiveInsertBody(body, node);
                }
            }
        }
        //Check bottom left
        else if (x <= midX && y > midY){
            if (root->bottomLeft == nullptr){
                coords boundTopLeft = {root->boundTopLeft.first, root->boundTopLeft.second+root->length/2};
                coords boundBottomRight = {root->boundBottomRight.first-root->length/2, root->boundBottomRight.second};
                root->bottomLeft = new Node(body.id, body.position, body.mass, root->length/2, boundTopLeft, boundBottomRight, true);
            }
            else{
                Node* node = root->bottomLeft;
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
                    updateMassAndCentreofMass(body, node);
                    recursiveInsertBody(body, node);
                }
            }
        }
        //Check bottom right
        else if (x > midX && y > midY){
            if (root->bottomRight == nullptr){
                coords boundTopLeft = {root->boundTopLeft.first+root->length/2, root->boundTopLeft.second+root->length/2};
                root->bottomRight = new Node(body.id, body.position, body.mass, root->length/2, boundTopLeft, root->boundBottomRight, true);
            }
            else{
                Node* node = root->bottomRight;
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
                    updateMassAndCentreofMass(body, node);
                    recursiveInsertBody(body, node);
                }
            }
        }
    }

    double calcDistance(coords a, coords b){ // distance = sqrt((x2-x1)^2+(y2-y1)^2)
        double distanceX {a.first-b.first};
        double distanceY {a.second-b.second};
        return std::sqrt(distanceX*distanceX + distanceY*distanceY);
    }

    std::pair<double, double> attraction(const Body& body, const Node* node){
        double distanceX {node->centreOfMass.first-body.position.first};
        double distanceY {node->centreOfMass.second-body.position.second};
        double distance {std::sqrt(distanceX*distanceX + distanceY*distanceY + EPSILON*EPSILON)};

        // Gravitational Force: F = GMm/r^2, G = Gravitional constant
        double force {G * body.mass * node->mass / (distance*distance)};
        // double angle {std::atan2(distanceY, distanceX)};
        // double forceX {std::cos(angle) * force};
        // double forceY {std::sin(angle) * force};
        double forceX {(distanceX / distance) * force};
        double forceY {(distanceY / distance) * force};

        return {forceX, forceY};
    }

    void recursiveForceCalculation(const Body& body, const Node* node, double& totalFx, double& totalFy){
        if (node->isLeaf){
            if (node->bodyId==body.id){
                return;
            }
            leafCount++;
            auto [fx, fy] = attraction(body, node);
            totalFx += fx;
            totalFy += fy;
        }
        // if less than theta use approximation as one body
        else if((node->length/calcDistance(body.position, node->centreOfMass)) < theta) {
            approxCount++;
            auto [fx, fy] = attraction(body, node);
            totalFx += fx;
            totalFy += fy;
        }
        else{
            recurseCount++;
            if (node->topLeft) recursiveForceCalculation(body, node->topLeft, totalFx, totalFy);
            if (node->topRight) recursiveForceCalculation(body, node->topRight, totalFx, totalFy);
            if (node->bottomLeft) recursiveForceCalculation(body, node->bottomLeft, totalFx, totalFy);
            if (node->bottomRight) recursiveForceCalculation(body, node->bottomRight, totalFx, totalFy);
        }
    }

    

public:
    size_t leafCount = 0;
    size_t approxCount = 0;
    size_t recurseCount = 0;

    //constructor
    QuadTree(const std::vector<Body>& bodies) : bodies(bodies), root(nullptr) {}

    ~QuadTree(){
        delete root;
    }

    void insertBody(const Body& body){
        if (root == nullptr){
            Node* newNode = new Node(body.id, body.position, body.mass, SIMSIZE, {0,0}, {SIMSIZE, SIMSIZE}, true);
            root = newNode;
            return;
        } else if(root->isLeaf){
            subdivide(body, root);
        } else{
            recursiveInsertBody(body, root);
            updateMassAndCentreofMass(body, root);
        }
    }

    void updateBodyPosition(Body& body){
        double totalFx {};
        double totalFy {};
        recursiveForceCalculation(body, root, totalFx, totalFy);

        body.xVel += totalFx / body.mass * TIMESTEP;
        body.yVel += totalFy / body.mass * TIMESTEP;

        body.position.first += body.xVel * TIMESTEP;
        body.position.second += body.yVel * TIMESTEP;
    }
};


int main()
{
	sf::RenderWindow window( sf::VideoMode( { WINDOWSIZE, WINDOWSIZE } ), "Barnes-Hut N Body Simulation" );

    sf::Clock clock;
    float fps;
    float averageFPS{0};
    float numFPSCounted{0};
    Timer t;

    std::mt19937 mt{std::random_device{}()};

    // --------- Change these - replace with other method of creating starting positions ---------
    // std::uniform_int_distribution createRandomPosition{0, static_cast<int>(SIMSIZE)};
    // std::uniform_real_distribution createRandomPosition{0.0, STARTINGSIZE}; // random position in starting sim
    // std::uniform_real_distribution createRandomRadius{0.0, STARTINGSIZE/2};
    std::uniform_real_distribution createRandomRadius{0.0, 1.0};
    std::uniform_real_distribution createRandomAngle{0.0, 2*PI};
    // std::uniform_real_distribution createRandomVelocity{-1*STARTINGVELOCITYRANGE, STARTINGVELOCITYRANGE};
    // --------------------------------

    sf::CircleShape bodyShape {};
    bodyShape.setRadius(2);
    bodyShape.setFillColor(COLOUR);
    bodyShape.setOrigin({2.0f,2.0f});

    std::vector<Body> bodies {};
    bodies.reserve(NUMOFBODIES);

    QuadTree* quadtree;
    quadtree = new QuadTree(bodies);

    double totalMass {STARTMASS*NUMOFBODIES};
    double offset {SIMSIZE/2};

    for (int i=0; i<NUMOFBODIES; i++){

        // --------- Change these - replace with other method of creating starting positions ---------
        // double radius {createRandomRadius(mt)};
        double radius {STARTINGSIZE/2 * std::sqrt(createRandomRadius(mt))};
        double angle {createRandomAngle(mt)};

        double x {radius*std::cos(angle)+offset};
        double y {radius*std::sin(angle)+offset};

        std::uniform_real_distribution createRandomVelocity{0.0, std::sqrt(2*G*totalMass/radius)};

        double xVel {0.0};
        double yVel {0.0};

        // double velocity {createRandomVelocity(mt)};
        // double direction {createRandomAngle(mt)};
        // double xVel {velocity*std::cos(direction)};
        // double yVel {velocity*std::sin(direction)};

        Body tempBody {i, {x, y}, STARTMASS, xVel, yVel};

        // Make bodies orbit centre
        // double dx = tempBody.position.first - SIMSIZE/2;
        // double dy = tempBody.position.second - SIMSIZE/2;
        // double dist = std::sqrt(dx*dx + dy*dy);
        // //speed required to maintain orbit, GM/r^2 = mv^2/r, v=sqrt(GM/R)
        // double fractionInterior = (dist * dist) / ((STARTINGSIZE/2) * (STARTINGSIZE/2));
        // double orbitalSpeed = std::sqrt(G * totalMass * fractionInterior / dist);

        // //perpindicular vector of (dx,dy) is (-dy,dx)
        // tempBody.xVel = -dy/dist * orbitalSpeed;
        // tempBody.yVel =  dx/dist * orbitalSpeed;
        // --------------------------------

        bodies.push_back(tempBody);
        quadtree->insertBody(bodies.back());
    }

    // window.setFramerateLimit(60); // turn on for constant frame rate
	while ( window.isOpen() )
	{
		while ( const std::optional event = window.pollEvent() )
		{
			if ( event->is<sf::Event::Closed>() ){
                std::cout << "Average FPS: " << averageFPS; //toggle
				window.close();
            }
		}

        // ------- toggle -------
        float currentTime = clock.restart().asSeconds();
        fps = 1.0f / currentTime;
        averageFPS = {(averageFPS*numFPSCounted+static_cast<int>(fps))/(numFPSCounted+1)};
        ++numFPSCounted;
        
        // std::cout << "FPS: " << static_cast<int>(fps) << "\n";

        // ----------------------

		window.clear();

        t.reset(); // toggle

        // std::cout << bodies.size() << '\n';
        // for (Body& body: bodies){
        //     if (!body.active) continue;

        //     if (body.position.first < 0 || body.position.second < 0 || 
        //         body.position.first > SIMSIZE || body.position.second > SIMSIZE){
        //             body.active = false;
        //             continue;
        //         }

        //     quadtree->updateBodyPosition(body);
        //     bodyShape.setPosition({static_cast<float>(body.position.first * SCALE), static_cast<float>(body.position.second * SCALE)});
        //     window.draw(bodyShape);
        // }

        // Not used as no discernable performace benefit, either due to gcc not implementing it or just equivalent performances
        std::for_each(std::execution::par_unseq, bodies.begin(), bodies.end(), 
                        [quadtree](Body& body){
                            if (!body.active) return;

                            if (body.position.first < 0 || body.position.second < 0 || 
                                body.position.first > SIMSIZE || body.position.second > SIMSIZE){
                                    body.active = false;
                                    return;
                                }
                            quadtree->updateBodyPosition(body);
                        });

        for (Body& body: bodies){
            if (!body.active) continue;
            bodyShape.setPosition({static_cast<float>(body.position.first * SCALE), static_cast<float>(body.position.second * SCALE)});
            window.draw(bodyShape);
        }

        window.display();

        std::cout << "Time to update bodies positions and draw to screen: " << t.elapsed() << "seconds\n"; //toggle
        t.reset();

        // std::cout << "leafCount: " << quadtree->leafCount << "\n";
        // std::cout << "approxCount: " << quadtree->approxCount << "\n";
        // std::cout << "recurseCount: " << quadtree->recurseCount << "\n";

        delete quadtree;
        std::cout << "Time to delete bodies: " << t.elapsed() << "seconds\n"; //toggle
        t.reset(); //toggle

        quadtree = new QuadTree(bodies);
        for (Body& body: bodies){
            if (!body.active) continue;
            quadtree->insertBody(body);
        }
        std::cout << "Time to create quadtree and insert bodies: " << t.elapsed() << "seconds\n"; //toggle
        t.reset();//toggle

        // remove body from bodies vector if outside range of sim
        // bodies.erase(std::remove_if(bodies.begin(), bodies.end(), [](const Body& body){
        //     return body.position.first < 0 || body.position.second < 0 
        //         || body.position.first > SIMSIZE || body.position.second > SIMSIZE;
        // }), bodies.end());
	}
}