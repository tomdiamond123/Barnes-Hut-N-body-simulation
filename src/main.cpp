#include <SFML/Graphics.hpp>
#include <iostream>
#include <random>

using coords = std::pair<double, double>;

double SIMSIZE {100*(3e8*365.25*24*60*60)};
int NUMOFBODIES {100};


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
};


class QuadTree {
private:
    //pointer to the root of the tree
    Node* root;

    void updateMassAndCentreofMass(const Body& body, Node* node){
        auto [nodeX, nodeY] = node->centreOfMass;
        node->mass += body.mass;
        nodeX = (nodeX*node->mass + body.position.first*body.mass)/node->mass;
        nodeY = (nodeY*node->mass + body.position.second*body.mass)/node->mass;
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
                updateMassAndCentreofMass(body, node);
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
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
                updateMassAndCentreofMass(body, node);
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
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
                updateMassAndCentreofMass(body, node);
                if (node->isLeaf){
                    subdivide(body, node);
                }
                else{
                    recursiveInsertBody(body, node);
                }
            }
        }
    }

public:
    //constructor
    QuadTree() : root(nullptr) {}

    void insertBody(const Body& body){
        

        if (root == nullptr){
            Node* newNode = new Node(body.position, body.mass, SIMSIZE, {0,0}, {SIMSIZE, SIMSIZE}, true);
            root = newNode;
            return;
        } else{
            recursiveInsertBody(body, root);
        }
    }


};


int main()
{
	QuadTree quadtree {};
    std::mt19937 mt{std::random_device{}()};
    // std::uniform_int_distribution createRandomPosition{0, static_cast<int>(SIMSIZE)};
    std::uniform_real_distribution createRandomPosition{0.0, SIMSIZE};

    for (int i; i<NUMOFBODIES; i++){
        Body tempBody {{createRandomPosition(mt), createRandomPosition(mt)}, 1e30, 0.0, 0.0};
        quadtree.insertBody(tempBody);
    }
	
	sf::RenderWindow window( sf::VideoMode( { 200, 200 } ), "SFML works!" );
	sf::CircleShape shape( 100.f );
	shape.setFillColor( sf::Color::Green );

	while ( window.isOpen() )
	{
		while ( const std::optional event = window.pollEvent() )
		{
			if ( event->is<sf::Event::Closed>() )
				window.close();
		}

		window.clear();
		window.draw( shape );
		window.display();
	}
}
