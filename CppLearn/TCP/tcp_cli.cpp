#include <cstdio>
#include <iostream>

#include "tcp_socket.hpp"

using namespace std;

int main(int argc, const char *argv[])
{
    if(argc != 3)
	{
		cout << "Usage:./tcp_cli ip port" << endl;
		return -1;
	}

    string srv_ip = argv[1];
    uint16_t src_port = stoi(argv[2]);

    TcpSocket tcpsocket;
    CHECK_RET(tcpsocket.Socket());

    CHECK_RET(tcpsocket.Connect(srv_ip, src_port));

    while(1)
    {
        string buf;
        cout << "Client Say: ";
        cin >> buf;
        tcpsocket.Send(buf);

        buf.clear();
        tcpsocket.Receive(&buf);
        cout << "Server Say: " << buf << endl;
    }
    
    tcpsocket.Close();
    return 0;
}
