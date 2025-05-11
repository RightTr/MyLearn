#include <cstdio>
#include <iostream>
#include <signal.h>
#include <sys/wait.h>

#include "tcp_socket.hpp"

using namespace std;

/*  处理子进程退出时产生的 SIGCHLD 信号
    防止僵尸进程（zombie process）*/
void sigcb(int no)
{
    while(waitpid(-1, NULL, WNOHANG) > 0);
}

int main(int argc, const char *argv[])
{
	if (argc != 3)
	{
		cout << "Usage:./tcp_srv ip port" << endl;
		return -1;
	}

    signal(SIGCHLD, sigcb); // 子进程退出处理方式
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

        int pid = fork(); //创建子进程
        /* 如果 pid == 0：说明这是新建的子进程。
           如果 pid > 0：还是原来的主进程 */
        if(pid == 0)
        {
            while(1)
            {
                string buf;
                new_tcpsocket.Receive(&buf);
                cout << "Client Say: " << buf << endl;

                buf.clear();
                cout << "Server Say: ";
                cin >> buf;
                new_tcpsocket.Send(buf);
            }
            new_tcpsocket.Close();
            exit(0); // 退出
        }
        new_tcpsocket.Close();
    }

    tcpsocket_list.Close();
    return 0;
}