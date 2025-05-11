#include <cstdio>
#include <iostream>

#include "tcp_socket.hpp"

int main(int argc, const char *argv[])
{
	if (argc != 3)
	{
		cout << "Usage:./tcp_srv ip port" << endl;
		return -1;
	}

    string ip = argv[1];
    uint16_t port = stoi(argv[2]);
    
    TcpSocket tcpsocket_list;
    CHECK_RET(tcpsocket_list.Socket());

    CHECK_RET(tcpsocket_list.Bind(&ip, port));

    CHECK_RET(tcpsocket_list.Listen());

    while(1)
    {
        TcpSocket new_tcpsocket;
        if(!tcpsocket_list.Accept(&new_tcpsocket))
        {
            continue; //服务端不能因为获取一个新建套接字失败就退出
        }
        string buf;
        new_tcpsocket.Receive(&buf);
        cout << "Client Say: " << buf << endl;

        buf.clear();
        cout << "Server Say: ";
        cin >> buf;
        new_tcpsocket.Send(buf);
    }

    tcpsocket_list.Close();
    return 0;
}

