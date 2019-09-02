/** \file magTCSsim.cpp
  * \author Jared R. Males
  * \brief Main program for magTCSsim
  */

//***********************************************************************//
// Copyright 2018 Jared R. Males (jaredmales@gmail.com)
//
// This file is part of magTCSsim.
//
// magTCSsim is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// magTCSsim is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with magTCSsim.  If not, see <http://www.gnu.org/licenses/>.
//***********************************************************************//


#include "magtelescope.hpp"

#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <pthread.h>


///Structure to pass to threads.
struct mtsfile
{
   magtelescope *mts;
   int fd;
   int port;
};

/* set up a listen socket
 */
int lsocket (int port)
{
   unsigned short sport;
   int s;
   struct  sockaddr_in sockaddr;
   int on;

   /* start listening for connection
    */
   sport = htons((short)port);
   if ((s = socket (AF_INET, SOCK_STREAM, 0)) < 0) 
   {
       perror ("socket");
       return (-1);
   }
   on = 1;
   if (setsockopt(s,SOL_SOCKET,SO_REUSEADDR,(char *)&on,sizeof(on)) < 0) 
   {
       perror ("setsockopt");
       return (-1);
   }
   memset ((char *)&sockaddr, 0, sizeof(sockaddr));
   sockaddr.sin_family  = AF_INET;
   sockaddr.sin_port    = sport;
   sockaddr.sin_addr.s_addr = INADDR_ANY;
   if (bind (s, (struct sockaddr *)&sockaddr, sizeof(sockaddr)) < 0) 
   {
       perror ("bind");
       return (-1);
   }
   if (listen (s, 10) < 0) 
   {
       perror ("listen");
       return (-1);
   }
   return (s);
}

void * listen_telescope_socket(void *vmf)
{
   mtsfile mtf;
   char string[MAGTCS_INPUT_SIZE];

   mtf.mts = ((mtsfile *) vmf)->mts;
   mtf.fd = ((mtsfile *) vmf)->fd;
   mtf.port = ((mtsfile*) vmf)->port;

   FILE *fp;
      
   if ((fp = fdopen (mtf.fd, "r+")) == NULL) 
   {
      exit (-6);
   }
      
   while (!tcs_fgets(string, MAGTCS_INPUT_SIZE, fp)) 
   {
      mtf.mts->process_string(std::string(string), fp, mtf.port == 5800 ); /// \todo status only port
   }
   
   std::cerr << "Closing\n";
   fclose(fp);
}

int network_accept_any( int & sockNo,
                        int socks[], 
                        unsigned int count,
                        struct sockaddr *addr, 
                        socklen_t *addrlen
                      ) 
{
   fd_set readfds;
   int maxfd, fd;
   unsigned int i;
   int status;

   FD_ZERO(&readfds);
   maxfd = -1;
   for(i = 0; i < count; i++) 
   {
      FD_SET(socks[i], &readfds);
      
      if (socks[i] > maxfd) maxfd = socks[i];
   }
    
   status = select(maxfd + 1, &readfds, NULL, NULL, NULL);
    
   if (status < 0)
   {
      return -1;
   }
    
   fd = -1;
   for (i = 0; i < count; i++)
   {
      if (FD_ISSET(socks[i], &readfds)) 
      {
         fd = socks[i];
         sockNo = i;
         break;
      }
   }  
   
   if (fd == -1)
      return -1;
   else
      return accept(fd, addr, addrlen);
}

int main (int argc, char *argv[])
{
   
   int s;
   
   magtelescope mts;

   mts.set_typval();

   mts.start_simulating();

   mts.start_tracking();

   
   int socks[2];
   int ports[2];
   
   ports[0] = 5800;
   ports[1] = 5811;
   
   if ((socks[0] = lsocket (ports[0])) < 0) 
   {
      fprintf (stderr, "lsocket (5800)\n");
      exit (-4);
   }
   
   if ((socks[1] = lsocket (ports[1])) < 0) 
   {
      fprintf (stderr, "lsocket (5811)\n");
      exit (-4);
   }
   
   while(1)
   {
      int fd;
      int sockNo;
      
      fd = network_accept_any(sockNo, socks, 2,(struct sockaddr *)0, (socklen_t *)0);
   
      mtsfile mtf;
       
      mtf.mts = &mts;
      mtf.fd = fd;
      mtf.port = ports[sockNo];
      
      std::cerr << "Connected: " << sockNo << " " << ports[sockNo] << " " << fd  << "\n";
      
      pthread_t thread_id;
      pthread_create(&thread_id, NULL, &listen_telescope_socket, (void *) &mtf);
   }
   
   exit (0);

}
