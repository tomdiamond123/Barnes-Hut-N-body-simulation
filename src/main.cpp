#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <algorithm>

using coords = std::pair<double, double>;

const double SIMSIZE {(3e8*365.25*24*60*60)};
const double SCALE {1000/SIMSIZE};
const int NUMOFBODIES {1000};
const double theta = 0.5;
const double G {6.67438e-11};
const double TIMESTEP = 3600*24*365;


struct Body{
    coords position{};
    double mass{};
    double xVel{};
    double yVel{};
};


class Node {
    public:
        coords centreOfMass{};
        double mass{};
        double length{};

        coords boundTopLeft{};
        coords boundBottomRight{};

        bool isLeaf{};

        Node* topLeft;
        Node* topRight;
        Node* bottomLeft;
        Node* bottomRight;

        // constructor
        Node(coords centreOfMass, double mass, double length, coords boundTopLeft, coords boundBottomRight, bool isLeaf) 
        : centreOfMass(centreOfMass), mass(mass), length(length), boundTopLeft(boundTopLeft), boundBottomRight(boundBottomRight), 
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
        Body oldBody {node->centreOfMass, node->mass, 0.0, 0.0};
        updateMassAndCentreofMass(body, node);
        recursiveInsertBody(oldBody, node);
        recursiveInsertBody(body, node);
    }

    void recursiveInsertBody(const Body& body, Node* root){
        auto [x, y] = body.position;

        double midX = (root->boundTopLeft.first + root->boundBottomRight.first) / 2.0;
        double midY = (root->boundTopLeft.second + root->boundBottomRight.second) / 2.0;

        //Check top left
        if (x <= midX && y <= midY){
            if (root->topLeft == nullptr){
                coords boundBottomRight = {root->boundBottomRight.first-root->length/2, root->boundBottomRight.second-root->length/2};
                root->topLeft = new Node(body.position, body.mass, root->length/2, root->boundTopLeft, boundBottomRight, true);
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
                root->topRight = new Node(body.position, body.mass, root->length/2, boundTopLeft, boundBottomRight, true);
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
                root->bottomLeft = new Node(body.position, body.mass, root->length/2, boundTopLeft, boundBottomRight, true);
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
                root->bottomRight = new Node(body.position, body.mass, root->length/2, boundTopLeft, root->boundBottomRight, true);
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

    double calcDistance(coords a, coords b){
        double distanceX {a.first-b.first};
        double distanceY {a.second-b.second};
        return std::sqrt(distanceX*distanceX + distanceY*distanceY);
    }

    std::pair<double, double> attraction(const Body& body, const Node* node){
        double distance {};
        double distanceX {node->centreOfMass.first-body.position.first};
        double distanceY {node->centreOfMass.second-body.position.second};
        distance = calcDistance(body.position, node->centreOfMass);

        double force {G * body.mass * node->mass / (distance*distance)};
        double angle {std::atan2(distanceY, distanceX)};
        double forceX {std::cos(angle) * force};
        double forceY {std::sin(angle) * force};

        return {forceX, forceY};
    }

    void recursiveForceCalculation(const Body& body, const Node* node, double& totalFx, double& totalFy){
        if (node->isLeaf){
            auto [fx, fy] = attraction(body, node);
            totalFx += fx;
            totalFy += fy;
        }
        else if((node->length/calcDistance(body.position, node->centreOfMass)) < theta) {
            auto [fx, fy] = attraction(body, node);
            totalFx += fx;
            totalFy += fy;
        }
        else{
            recursiveForceCalculation(body, node->topLeft, totalFx, totalFy);
            recursiveForceCalculation(body, node->topRight, totalFx, totalFy);
            recursiveForceCalculation(body, node->bottomLeft, totalFx, totalFy);
            recursiveForceCalculation(body, node->bottomRight, totalFx, totalFy);
        }
    }

    

public:
    //constructor
    QuadTree() : root(nullptr) {}

    ~QuadTree(){
        delete root;
    }

    void insertBody(const Body& body){
        if (root == nullptr){
            Node* newNode = new Node(body.position, body.mass, SIMSIZE, {0,0}, {SIMSIZE, SIMSIZE}, true);
            root = newNode;
            return;
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
	sf::RenderWindow window( sf::VideoMode( { 1000, 1000 } ), "Barnes-Hut N Body Simulation" );

    QuadTree* quadtree;
    quadtree = new QuadTree{};
    std::mt19937 mt{std::random_device{}()};
    // std::uniform_int_distribution createRandomPosition{0, static_cast<int>(SIMSIZE)};
    std::uniform_real_distribution createRandomPosition{0.0, SIMSIZE};

    sf::CircleShape bodyShape {};
    bodyShape.setRadius(2);
    bodyShape.setFillColor(sf::Color::Yellow);
    bodyShape.setOrigin({2.0f,2.0f});

    std::vector<Body> bodies {};

    for (int i; i<NUMOFBODIES; i++){
        Body tempBody {{createRandomPosition(mt), createRandomPosition(mt)}, 1e30, 0.0, 0.0};
        quadtree->insertBody(tempBody);
        bodies.push_back(tempBody);


        // bodyShape.setPosition({static_cast<float>(tempBody.position.first * SCALE), static_cast<float>(tempBody.position.second * SCALE)});
        // window.draw(bodyShape);
    }

    // window.setFramerateLimit(60);
	while ( window.isOpen() )
	{
		while ( const std::optional event = window.pollEvent() )
		{
			if ( event->is<sf::Event::Closed>() )
				window.close();
		}

		window.clear();

        for (Body& body: bodies){
            quadtree->updateBodyPosition(body);
            bodyShape.setPosition({static_cast<float>(body.position.first * SCALE), static_cast<float>(body.position.second * SCALE)});
            window.draw(bodyShape);
        }
        window.display();


        delete quadtree;
        quadtree = new QuadTree{};
        for (Body& body: bodies){
            quadtree->insertBody(body);
        }

        bodies.erase(std::remove_if(bodies.begin(), bodies.end(), [](const Body& body){
            return body.position.first < 0 || body.position.second < 0 
                || body.position.first > SIMSIZE || body.position.second > SIMSIZE;
        }), bodies.end());
	}
}