#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct nodeValue {
    int node_id;
    int n_dims;
    double radius;

    double * center_coordinates;
};

struct node {
  struct nodeValue item;
  struct node* left;
  struct node* right;
};

// Postorder traversal
void postorderTraversal(struct node* root) {
  if (root == NULL) return;
  postorderTraversal(root->left);
  postorderTraversal(root->right);

  int j;
  int left_child_id = root->left == NULL ? -1 : root->left->item.node_id;
  int right_child_id = root->right == NULL ? -1 : root->right->item.node_id;;
  printf("%d %d %d %f", root->item.node_id, left_child_id, right_child_id, root->item.radius);

  for(j = 0; j < root->item.n_dims && root->item.center_coordinates != NULL; j++){
      printf(" %f", root->item.center_coordinates[j]);
  }


  printf("\n");
}

// Postorder traversal
void clean(struct node* root) {
    if (root == NULL) return;
    clean(root->left);
    clean(root->right);

    free(root->item.center_coordinates);
    free(root);
}

// Create a new Node
struct node* createNode(struct nodeValue value) {
  struct node* newNode = malloc(sizeof(struct node));
  newNode->item = value;
  newNode->left = NULL;
  newNode->right = NULL;

  return newNode;
}