/*
*  Copyright (C) 2007 Bernd Machenschalk, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*
  Author Cokelaer Thomas
  Functions to manipulate linked list.
*/


#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdio.h>
#include <lal/FindRoot.h>

/* for debugging only
void
print_list(CellList *head)
{
  if (head == NULL){
    printf("\n");
  }
  else {
    printf(" %d", head->id);
    print_list(head->next);
  }
}
*/

UINT4 LALListLength(CellList *list)
{
  UINT4 count = 0;
  while (list != NULL){
    count++;
    list = list->next;
  }
  return count;
}



void LALListAppend(
	CellList **headRef,
	INT4 id)
{
  CellList *current;

  if ((current = malloc(sizeof(*current))) == NULL) {
    {
      printf("Error with malloc\n");
      exit(0);
    }
  }
  current->id = id;
  current->next = *headRef;
  *headRef = current;
}



/*Is it used somewhere ? delete everything in principle*/
/*
void DeleteList(CellList **headRef)
{
  CellList *tmp;

  while (headRef !=NULL){
    tmp = (*headRef)->next;
    free(headRef);
    (*headRef) = tmp;
  }

}
*/


void LALListDelete(CellList **headRef, INT4 id)
{
  CellList *ptr  = NULL;
  CellList *prev = NULL;


  for (ptr = *headRef; ptr != NULL; ptr= ptr->next){
    if (id == ptr->id)
      break;

    prev = ptr;
  }

  if (ptr==NULL)
    return ;

  if (prev!=NULL){
    prev->next = ptr->next;
  }
  else{
    *headRef = ptr->next;
  }

  /* free the data here if needed such as free(ptr->id); */
  free(ptr);

  return ;



}
