#include <time.h>
#include <iostream>
#include <queue>
#include <string>
#include <algorithm>
using namespace std;
#define M 2//m*2=key����m*2+1=ֵ��//��������࣬ɾ��һ���ڵ㲻������(M+1)/2,��zhiΪ5����������2
string key_words[1001];//5����,���4���ؼ���
class Inter_Node;

class Node
{
public:
    Node();
    virtual ~Node();
    Node* GetBrother(int& flag);
    Inter_Node* Parent;
    int key[M * 2];//4
    int count;
    int isLeaf;
    void Print();
};
Node::~Node() {}

class Inter_Node : public Node
{
public:
    Inter_Node();
    virtual ~Inter_Node();//�麯����������
    bool Insert(int value, Node* pNode);
    bool Delete(int value);
    int Split(Inter_Node* pNode, int key);//����
    bool Merge(Inter_Node* pNode);
    bool Slib(Inter_Node* pNode);
    Node* Child[M * 2 + 1];//5������·��
};
Inter_Node::~Inter_Node()
{
    for (int i = 0; i < M * 2 + 1; i++)
        Child[i] = NULL;
}

//Ҷ�ӽ��
class Leaf_Node : public Node
{
public:
    Leaf_Node();
    virtual ~Leaf_Node();
    bool Insert(int value);
    bool Delete(int value);
    int Split(Leaf_Node* pNode);
    bool Merge(Leaf_Node* pNode);
    Leaf_Node* Pre_Node;
    Leaf_Node* Next_Node;
};
Leaf_Node::~Leaf_Node() {}

//B+��
class Bplus
{
public:
    Bplus();
    virtual ~Bplus();
    bool Search(int data, string& sPath);
    bool Insert(int data);
    bool Delete(int data);
    void Print();

protected:
    Leaf_Node* Find(int data);
    bool Add_Node(Inter_Node* pNode, int key, Node* New_Node);
    bool Remove_Node(Inter_Node* pNode, int key);
    Node* Root;
};
Bplus::~Bplus() {}

//��㴴������
Node::Node()
{
    isLeaf = true;
    count = 0;
    Parent = NULL;
}
//Ҷ�ӽ�㴴������
Leaf_Node::Leaf_Node()
{
    isLeaf = true;
    Pre_Node = NULL;
    Next_Node = NULL;
}
//�м��㴴������
Inter_Node::Inter_Node()
{
    isLeaf = false;
    for (int i = 0; i < M * 2 + 1; i++)
        Child[i] = NULL;
}
//Bplus��������
Bplus::Bplus()
{
    Root = NULL;
}

//�������ֵܽ��
Node* Node::GetBrother(int& flag)
{
    if (NULL == Parent)
        return NULL;

    Node* p = NULL;
    for (int i = 0; i <= Parent->count; i++)
    {
        if (Parent->Child[i] == this)
        {
            if (i == Parent->count)
            {
                p = Parent->Child[i - 1];//���ֵ� flag=1
                flag = 1;
            }
            else
            {
                p = Parent->Child[i + 1];//���ֵ� flag=2
                flag = 2;
            }
        }
    }
    return p;
}

//������
void Node::Print()
{
    for (int i = 0; i < count; i++)
    {
        cout << "(" << key[i] << "," << key_words[key[i]] << ")  ";
        if (i >= count - 1)
            cout << " | ";
    }
}

//Ҷ�ӽ��ķ���
int Leaf_Node::Split(Leaf_Node* p)
{
    int j = 0;
    for (int i = M; i < M * 2; i++, j++)//��ֵcopy���½ڵ�
        p->key[j] = this->key[i];//thisΪold node
    this->count = this->count - j;
    p->count = j;
    return p->key[0];
}

//Ҷ�ӽ��ɾ��
bool Leaf_Node::Delete(int value)
{
    bool found = false;
    int i = 0;
    for (; i < count; i++)
    {
        if (value == key[i])
        {
            found = true;
            break;
        }
    }
    if (false == found)
        return false;
    int j = i;
    for (; j < count - 1; j++)
        key[j] = key[j + 1];
    key[j] = 0;
    count--;
    return true;
}

//Ҷ�ӽ��Ĳ���
bool Leaf_Node::Insert(int value)
{
    int i = 0;
    for (; (value > key[i]) && (i < count); i++)//��˳��
    {
    }
    for (int j = count; j > i; j--)//�ƶ����ҵ�Ӧ�ò���Ĺؼ���λ��
        key[j] = key[j - 1];
    key[i] = value;//����ؼ���
    count++;
    return true;
}

//Ҷ������
Leaf_Node* Bplus::Find(int data)//dataΪ�ؼ���
{
    int i = 0;
    Node* p = Root; //?????bplus�ĸ�
    Inter_Node* q;  //?��???
    while (NULL != p)
    {
        if (p->isLeaf) //????????????
            break;
        for (i = 0; i < p->count; i++) //??????p??key?????p����Ҷ�ӣ�ѭ����count�Ľڵ�
        {
            if (data < p->key[i])
                break;
        }
        q = (Inter_Node*)p;
        p = q->Child[i];
    }

    return (Leaf_Node*)p;//�Ѹ�return,�����Ϊ��,��һ�����뺯�����ɵĽڵ㼴Ϊ��
}

//Ҷ�ӽ��,�ϲ�Ҷ�ӽ��
bool Leaf_Node::Merge(Leaf_Node* p)
{
    if (this->count + p->count > M * 2)//�������һ�����˵������Ҫ�ϲ�
        return false;
    for (int i = 0; i < p->count; i++)//����oldnode�Ĺؼ��ֶ����뵽bro��
        this->Insert(p->key[i]);
    return true;
}

//�м���Merge�ϲ�
bool Inter_Node::Merge(Inter_Node* p)
{
    key[count] = p->Child[0]->key[0];
    count++;
    Child[count] = p->Child[0];
    for (int i = 0; i < p->count; i++)
    {
        key[count] = p->key[i];
        count++;
        Child[count] = p->Child[i + 1];
    }
    return true;
}

//�м������
bool Inter_Node::Insert(int value, Node* New)
{
    int i = 0;
    for (; (i < count) && (value > key[i]); i++)//iָ��keyҪ�����λ��
    {
    }
    for (int j = count; j > i; j--)//Ų�����ط�
        key[j] = key[j - 1];
    for (int j = count + 1; j > i + 1; j--)//����keyֵ�ı䣬�����ƶ���
        Child[j] = Child[j - 1];
    key[i] = value;//�ؼ��ִ������׽ڵ�
    Child[i + 1] = New;//newnode�ŵ��÷ŵ�λ��
    New->Parent = this;
    count++;
    return true;
}

//�м������
int Inter_Node::Split(Inter_Node* p, int k)
{
    int i = 0, j = 0;
    if ((k > this->key[M - 1]) && (k < this->key[M]))//���ѵĵط����м�
    {
        for (i = M; i < M * 2; i++, j++)
            p->key[j] = this->key[i];//��������ֵ��brother
        j = 1;
        for (i = M + 1; i <= M * 2; i++, j++)
        {
            this->Child[i]->Parent = p;//���Ӹ��������ƶ�
            p->Child[j] = this->Child[i];
        }
        this->count = M;//�ؼ���������λһ��
        p->count = M;
        return k;
    }
    int pos = k < this->key[M - 1] ? (M - 1) : M;//��k��С���м�-1�Ƚϣ���λ��ǰ�滹���ں���ڵ�
    k = this->key[pos];//posΪ���ѵ�,��λΪǰ���Ǻ���ѵ�,���϶�Ϊ�м�ֵ
    j = 0;
    for (i = pos + 1; i < M * 2; i++, j++)//ǰ�ڵ㿼��ڵ�,�Ӳ����λ�÷֣������Ժ�ķŵ��½ڵ�
        p->key[j] = this->key[i];
    j = 0;
    for (i = pos + 1; i <= M * 2; i++, j++)//�������͸��ֵ�
    {
        this->Child[i]->Parent = p;
        p->Child[j] = this->Child[i];
    }
    this->count = pos;
    p->count = M * 2 - pos - 1;
    return k;
}

//�м���ɾ��
bool Inter_Node::Delete(int k)
{
    int i = 0;
    for (; (k >= key[i]) && (i < count); i++)
    {
    }
    for (int j = i - 1; j < count - 1; j++)
        key[j] = key[j + 1];
    k = i;
    for (; k < count; k++)
    {
        Child[k] = Child[k + 1];
    }
    Child[k] = NULL;
    count--;
    return true;
}

//�м���
bool Inter_Node::Slib(Inter_Node* p)
{
    int i, j;
    if (p->key[0] < this->key[0])
    {
        for (i = count; i > 0; i--)
            key[i] = key[i - 1];
        for (j = count + 1; j > 0; j--)
            Child[j] = Child[j - 1];
        key[0] = Child[0]->key[0];
        Child[0] = p->Child[p->count];
    }
    else
    {
        key[count] = p->Child[0]->key[0];
        Child[count + 1] = p->Child[0];
        for (i = 1; i < p->count - 1; i++)
            p->key[i - 1] = p->key[i];
        for (j = 0; j < p->count - 1; j++)
            p->Child[j] = p->Child[j + 1];
    }
    this->count++;
    p->count--;
    return true;
}

//B+����ӽ��
bool Bplus::Add_Node(Inter_Node* p, int k, Node* New_Node)
{
    if (NULL == p || p->isLeaf)
        return false;
    if (p->count < M * 2)//���ײ���
        return p->Insert(k, New_Node);
    Inter_Node* Brother = new Inter_Node;
    //Ҷ�ӽڵ��������ڵ�Ҳ���������
    int NewKey = p->Split(Brother, k);//NewKeyΪ��Ҫ��ȡ�����뵽root�ڵ��ֵ

    //ȷ����Ҫ����Ĺؼ��֣��ǲ��뵽���ѽڵ���ĸ�λ��
    if (p->count < Brother->count)
    {
        p->Insert(k, New_Node);
    }
    else if (p->count > Brother->count)
    {
        Brother->Insert(k, New_Node);
    }
    else
    {
        Brother->Child[0] = New_Node;
        New_Node->Parent = Brother;
    }
    Inter_Node* parent = (Inter_Node*)(p->Parent);
    if (NULL == parent)
    {
        parent = new Inter_Node();
        parent->Child[0] = p;
        parent->key[0] = NewKey;//newkeyΪ���Ѵ��أ�Ϊ������м�ֵ
        parent->Child[1] = Brother;
        p->Parent = parent;
        Brother->Parent = parent;
        parent->count = 1;
        Root = parent;
        return true;
    }
    return Add_Node(parent, NewKey, Brother);
}

//B+������data
bool Bplus::Search(int data, string& sPath)
{
    int i = 0;
    sPath = "����·��Ϊ: ";
    Node* p = Root;
    if (NULL == p)
        return false;
    Inter_Node* q;
    while (NULL != p)
    {
        if (p->isLeaf)
            break;
        for (i = 0; (i < p->count) && (data >= p->key[i]); i++)
        {
        }
        int k = i > 0 ? i - 1 : i;
        sPath += to_string(p->key[k]);
        sPath += "-->";
        q = (Inter_Node*)p;
        p = q->Child[i];
    }
    if (NULL == p)
        return false;
    sPath += to_string(p->key[0]);
    bool found = false;
    for (i = 0; i < p->count; i++)
    {
        if (data == p->key[i])
            found = true;
        //sPath += to_string(p->key[i]);
        //sPath += "-->";
    }
    if (found)
    {
        sPath += " OK";
    }
    else
    {
        sPath += " FAIL";
    }
    return found;
}

//B+���Ĳ���
bool Bplus::Insert(int data) //dataΪ����Ĺؼ���
{

    string a;
    if (true == Search(data, a))//����Ҫ�����ֵ
        return false;

    Leaf_Node* Old_Node = Find(data);//�ҵ���Ҫ�����Ҷ�ӽڵ㶨��Ϊoldnode

    if (NULL == Old_Node)
    {
        Old_Node = new Leaf_Node;//��Ϊ��
        Root = Old_Node;
    }

    if (Old_Node->count < M * 2) {//�пռ���룬ֱ�Ӳ��ȥ������
        return Old_Node->Insert(data);

    }

    Leaf_Node* New_Node = new Leaf_Node;//��������

    int k = Old_Node->Split(New_Node);//kΪ�½ڵ��һ���ؼ���

    Leaf_Node* OldNext = Old_Node->Next_Node;
    Old_Node->Next_Node = New_Node;//����Ҷ�ӽڵ�����
    New_Node->Next_Node = OldNext;
    New_Node->Pre_Node = Old_Node;

    if (NULL != OldNext)
        OldNext->Pre_Node = New_Node;

    if (data < k)//С��newnode key[0]����ǰ�棬��������
    {
        Old_Node->Insert(data);
    }
    else
    {
        New_Node->Insert(data);
    }
    Inter_Node* parent = (Inter_Node*)(Old_Node->Parent);

    if (NULL == parent)//��ʼ��parent����û�и���㣬�½�һ��
    {
        Inter_Node* New_Root = new Inter_Node;
        New_Root->Child[0] = Old_Node;
        New_Root->key[0] = k;
        New_Root->Child[1] = New_Node;
        Old_Node->Parent = New_Root;
        New_Node->Parent = New_Root;
        New_Root->count = 1;
        Root = New_Root;
        return true;
    }

    return Add_Node(parent, k, New_Node);//�������ֵ���߷��Ѹ��׽����µĽڵ�
}

//B+����ɾ��
bool Bplus::Delete(int data)
{
    Leaf_Node* Old_Node = Find(data); //��������
    if (NULL == Old_Node)//��Ϊ��
        return false;
    if (false == Old_Node->Delete(data)) //ɾ��
        return false;
    Inter_Node* parent = (Inter_Node*)(Old_Node->Parent);
    if (NULL == parent)
    {
        if (0 == Old_Node->count)//��������ɾ����û����ûkey
        {
            delete Old_Node;
            Root = NULL;
        }
        return true;
    }
    if (Old_Node->count >= M)
    {
        for (int i = 0; (i < parent->count) && (data >= parent->key[i]); i++)
        {
            if (parent->key[i] == data)//���Ҫɾ����key�ȸ�����������ֵҪ���ֱ��ɾ���������ȣ��͸�����һ��������
                parent->key[i] = Old_Node->key[0];
        }
        return true;
    }
    //��������Ҫ�ϲ����ֵ
    int flag = 1;
    Leaf_Node* Brother = (Leaf_Node*)(Old_Node->GetBrother(flag));
    int NewData = 0;
    if (Brother->count > M)//��ֵ
    {
        if (1 == flag)//���ֵ�
        {
            NewData = Brother->key[Brother->count - 1];//Ҫ�����ߵ�����
        }
        else//���ֵ�
        {
            NewData = Brother->key[0];
        }
        Old_Node->Insert(NewData);
        Brother->Delete(NewData);
        //�滻parent�е�keyֵ
        if (1 == flag)
        {
            for (int i = 0; i <= parent->count; i++)//�����ֵܽ�ֵ
            {
                if (parent->Child[i] == Old_Node && i > 0)
                    parent->key[i - 1] = Old_Node->key[0];
            }
        }
        else
        {
            for (int i = 0; i <= parent->count; i++)//�����ֵܽ�ֵ
            {
                if (parent->Child[i] == Old_Node && i > 0)
                    parent->key[i - 1] = Old_Node->key[0];
                if (parent->Child[i] == Brother && i > 1)
                    parent->key[i - 1] = Brother->key[0];
            }
        }
        return true;
    }
    int NewKey = 0;
    if (1 == flag)//�޷���ֵ���ϲ�
    {
        Brother->Merge(Old_Node);
        NewKey = Old_Node->key[0];//���Ҫɾ���ĸ������key
        Leaf_Node* OldNext = Old_Node->Next_Node;//��������ֵ�
        Brother->Next_Node = OldNext;
        if (NULL != OldNext)
            OldNext->Pre_Node = Brother;
        delete Old_Node;
    }
    else
    {
        Old_Node->Merge(Brother);
        NewKey = Brother->key[0];
        Leaf_Node* OldNext = Brother->Next_Node;
        Old_Node->Next_Node = OldNext;
        if (NULL != OldNext)
            OldNext->Pre_Node = Old_Node;
        delete Brother;
    }
    return Remove_Node(parent, NewKey);//�Ƴ�parent�����Ƴ�parent�йؼ��֣�
}

//Bplus �Ƴ����
bool Bplus::Remove_Node(Inter_Node* p, int k)
{
    if (false == p->Delete(k))
    {
        return false;
    }
    Inter_Node* parent = (Inter_Node*)(p->Parent);
    if (NULL == parent)
    {
        if (0 == p->count)
        {
            Root = p->Child[0];
            delete p;
        }
        return true;
    }
    if (p->count >= M)//���ײ��ϲ�
    {
        //ɾ��parent�еĹؼ���
        for (int i = 0; (i < parent->count) && (k >= parent->key[i]); i++)
        {
            if (parent->key[i] == k)//�����׵�parent����û��Ҫɾ���Ĺؼ��֣��о͸�������
                parent->key[i] = p->key[0];
        }
        return true;
    }
    //���׺ϲ�
    int flag = 1;
    Inter_Node* Brother = (Inter_Node*)(p->GetBrother(flag));
    if (Brother->count > M)//���׽�ֵ
    {
        p->Slib(Brother);
        if (1 == flag)
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == p && i > 0)
                    parent->key[i - 1] = p->key[0];
            }
        }
        else
        {
            for (int i = 0; i <= parent->count; i++)
            {
                if (parent->Child[i] == p && i > 0)
                    parent->key[i - 1] = p->key[0];
                if (parent->Child[i] == Brother && i > 0)
                    parent->key[i - 1] = Brother->key[0];
            }
        }
        return true;
    }
    //�ֵܽ�ֵ
    int NewKey = 0;
    if (1 == flag)
    {
        Brother->Merge(p);
        NewKey = p->key[0];
        delete p;
    }
    else
    {
        p->Merge(Brother);
        NewKey = Brother->key[0];
        delete Brother;
    }
    return Remove_Node(parent, NewKey);
}

//Bplus���
void Bplus::Print()
{
    Node* p = Root;
    if (NULL == p)
        return;
    Inter_Node* a;
    int H = 0;
    queue<Node*> q;
    queue<int> h;
    q.push(p);
    h.push(1);
    while (!q.empty())
    {
        p = q.front();
        if (H != h.front())
        {
            cout << endl;
            cout << H << endl;
            H = h.front();
        }
        q.pop();
        h.pop();
        p->Print();
        if (NULL != p && !p->isLeaf)
        {
            a = (Inter_Node*)p;
            for (int i = 0; i <= p->count; i++)
            {
                q.push(a->Child[i]);
                h.push(H + 1);
            }
        }
    }
}

//���ԣ�test1 ����B+tree
void test1(Bplus* bplus, int count)//countΪ��Ҫ�����ڵ㣬����ֵ��һ����
{
    srand((unsigned)time(NULL));
    //�������
    cout << " �������ݣ� ";
    for (int i = 1; i <= count; i++)//
    {
        //int x = i;
        int x = rand() % 1000 + 1;//�������x
        key_words[x] = " ";
        for (int i = 0; i < 2; i++)
        {
            int y = rand() % 26 + 65;//����һ�������
            int z = rand() % 26 + 97;
            key_words[x] += char(y);
            key_words[x] += char(z);
        }//�����ĸ��ַ�������
        cout << "[" << x << " "
            << "," << key_words[x] << "]" << " ";
        bplus->Insert(x);//x=key
    }
    cout << "B+������suceess!" << endl;
    cout << endl;
}

//���ԣ�test2��ѯ
void test2(Bplus* bplus, int data)
{
    string sPath;
    bplus->Search(data, sPath);
    cout << sPath;
    cout << endl;
}

//���ԣ�test3����
void test3(Bplus* bplus, int data, string s)
{
    bool success = bplus->Insert(data);
    key_words[data] = s;
    if (true == success)
    {
        cout << "OK" << endl;
    }
    else
    {
        cout << "FAIL" << endl;
    }
    cout << endl;
}

//���ԣ�test4ɾ��
void test4(Bplus* bplus, int data)
{
    bool success = bplus->Delete(data);
    if (true == success)
    {
        cout << "�ɹ�" << endl;
    }
    else
    {
        cout << "����" << endl;
    }
    cout << endl;
}

//���ԣ�test5��ӡ
void test5(Bplus* bplus)
{
    bplus->Print();
    cout << endl;
}

//���ԣ�test6�޸�
void test6(Bplus* bplus, int i, string j, int k, string l)
{
    if (i == k)
    {
        key_words[i] = l;
    }
    else
    {
        //i ��= j
        //ɾ��(i , j );
        bplus->Delete(i);
        //����(k , l )
        bplus->Insert(k);
        //��ʼ���������
        key_words[k] = l;
    }
    cout << endl;
}

//B+��������
int BPlusPlus()
{
    Bplus* bplus = new Bplus();
    int x = 1;
    int y = 0;
    string s = " ";
    int i, k;
    string j = " ";
    string l = " ";
    while (0 != x)
    {
        cout << "1.������Ҫ���ٸ��ڵ��B+tree" << endl;
        cout << "2.��Ҫ��ѯ " << endl;
        cout << "3.��Ҫ�������� " << endl;
        cout << "4.��Ҫɾ������" << endl;
        cout << "5.��Ҫ��ʾB+����״" << endl;
        cout << "6.��Ҫ�޸�����" << endl;
        cout << " \n";
        cin >> x;
        switch (x)
        {
        case 1:
            cout << "��������Ҫ���ٸ���㣺";
            cin >> y;
            test1(bplus, y);
            break;
        case 2:
            cout << "������Ҫ��ѯ�Ĺؼ�ֵkey A��";
            cin >> y;
            test2(bplus, y);
            break;
        case 3:
            cout << "������Ҫ����Ĺؼ�ֵkey A��B��";
            cin >> y >> s;
            test3(bplus, y, s);
            break;
        case 4:
            cout << "������Ҫɾ���Ĺؼ�ֵkey A��";
            cin >> y;
            test4(bplus, y);
            break;
        case 5:
            test5(bplus);
            break;
        case 6:
            cout << "������Ҫɾ����A��B��";
            cin >> i >> j;
            cout << "������Ҫ�����A��B��";
            cin >> k >> l;
            test6(bplus, i, j, k, l);
            break;
        default:
            delete bplus;
            return 0;
            break;
        }
    }
    return 0;
}

//��������������������������������
//��Ȩ����������ΪCSDN�����������ϲ�����ԭ�����£���ѭCC 4.0 BY - SA��ȨЭ�飬ת���븽��ԭ�ĳ������Ӽ���������
//ԭ�����ӣ�https://blog.csdn.net/weixin_51309915/article/details/122592023