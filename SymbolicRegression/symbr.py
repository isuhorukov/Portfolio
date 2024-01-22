import numpy as np
from sympy import simplify, cos, sin, exp, log, re, symbols, I, pi, Pow
from sympy import *

class node:
    def __init__(self, operation, left = None, right = None):
        self.operation  = operation
        self.par        = None
        self.dir        = ""
        self.left       = left
        self.right      = right

    def set_par(self, par):
        self.par = par

class symb_exp:

    def __init__(self, x_length = 3, number_of_mutations = 1, max_length_tree = 7):
        """
        :param x_length:            Размерность функции
        :param number_of_mutations: Количество мутаций за генерацию
        """
        pass
        self.x_length       = x_length
        self.NOM            = number_of_mutations
        self.MLT            = max_length_tree
        self.tree_length    = 0
        self.score          = -1
        self.score_MMS      = -1

        self.oper_list = []
        self.symb_list = []

        self.operators    = ["+", "*", "apow"]
        self.functions  = ["alog", "sin", "exp", "(-1.0)*", "cos"]

        self.oper_priority = {\
            "+" : 1,
            "*" : 2,
            "^re" : 3,
            "apow" : 3,
            "(-1.0)*" : 0
            }
        
        self.expression = []
        self.tree = []
        self.root = ""

        self.mut_string = ""
        self.eval_string = ""


    def get_symb_str(self, expression_get: str):
        """
        expression_get: Строка выражения в обычной записи
        """
        expr_list = expression_get.strip().split(" ")           # Изначально предполагается, что выражение записано через пробелы
        ind = 0
        self.expression = []
        while ind < len(expr_list):                             # Пока не дойдем до последнего выражения
            expr = expr_list[ind]                               # Получаем текущее выражение
            if (not(expr in self.operators or expr in self.functions) and expr != "(" and expr != ")" and expr != ","):
                self.expression.append(expr)                    #  /\ Здесь написано следующее: 
                                                                #  || Если выражение является символом (число и символ отождествляем),
                                                                #  || то добавим его в список польской записи

            if (expr in self.functions):
                self.oper_list.append(expr)

            if (expr in self.operators):                                    # Если выражение оказалось оператором
                if (len(self.oper_list) > 0):                               # и в списке операторов есть что-то
                    if (self.oper_list[-1] != "("):                         # и это не левая скобка
                        p1 = self.oper_priority[self.oper_list[-1]]         # и приоритет верхней операции в списке выше текущей
                        p2 = self.oper_priority[expr]                       # добавим операцию из списка в польскую запись
                        if(expr != "(" and p1 > p2):                        # <- Сравнение приоритетов
                            self.expression.append(self.oper_list.pop())    # <- Добавили в польскую запись


                self.oper_list.append(expr)                                 # Добавим текущее выражение в операции
            
            if (expr == ","):                                           # Если выражение это запятая (например для функции pow,max,min и пр.)
                expr_buff = self.oper_list[-1]                          # и пока верхняя операция в списке не равна левой скобке
                while(expr_buff != "("):                                # записываем все в польскую запись
                    self.expression.append(self.oper_list.pop())        # <- Достаем из списка

                    if(len(self.oper_list) != 0):
                        expr_buff = self.oper_list[-1]

            if (expr == "("):
                self.oper_list.append("(")


            if (expr == ")"):                                       # Если выражение правая скобка,
                expr_buff = self.oper_list[-1]                      # то пишем все из списка операций в польскую запись
                while(expr_buff != "("):                            # попутно проверяя, что верхний элемент есть
                    self.expression.append(self.oper_list.pop())    # Саму скобку надо вежливо убрать никуда не записывая
                    if(len(self.oper_list) != 0):                   #
                        expr_buff = self.oper_list[-1]              #
                self.oper_list.pop()                                # <- Убираем скобку


                if(len(self.oper_list) > 0):                                # Может оказаться, что скобки стояли 
                    if(self.oper_list[-1] in self.functions):               # для определения аргументов функции
                        self.expression.append(self.oper_list.pop())        # поэтому проверяем верхний элемент на принадлежность функции

            ind += 1


        while len(self.oper_list) > 0:                          # В конце алгоритма получается так, что остались операции без скобок и прочего
            if(self.oper_list[-1] == "("):                      # поэтому все операции надо записать в польскую запись,
                self.oper_list.pop()                            # а левую скобку удалить из списка
            else:                                               #
                self.expression.append(self.oper_list.pop())    #

        self.get_tree()
       
    def get_tree(self):
        # Функция получения дерева выражений из польской записи
        
        self.tree = []
        for expr in self.expression:                            # Идем по записи слева-направо          
            if(expr in self.operators):                         # Если элемент операция, то достаем два вверхних элемента, к которым применяется операция
                                                                # 
                right = self.tree.pop()                         # Достаем два вверхних элемента
                left = self.tree.pop()                          # Эти два элемента связаны операцией
                                                                # 
                self.tree.append(node(expr, left, right))       # Добавляем в дерево узел, с операцией expr и двумя узлами left и right
                self.tree[-1].left.set_par(self.tree[-1])
                self.tree[-1].right.set_par(self.tree[-1])
            else:
                if(expr in self.functions):
                                                                # Если элемент функция, то достаем только один элемент
                    right = self.tree.pop()
                                                                # Добавляем в дерево узел, с функцией expr и двумя узлами None и right
                    self.tree.append(node(expr,None, right))
                    self.tree[-1].right.set_par(self.tree[-1])
                else:
                                                                # Если элемент не функция и не операция, значит это число, либо символ (здесь и далее говорим просто символ)
                    self.tree.append(node(expr))
        

        self.root = self.tree[-1]                   # Вверхний элемент списка, нам достаточно хранить только корень дерева
        self.tree_length = self.get_tree_length()   # Сразу получаем длину дерева

                                    # Получаем две строки,
                                    # Первая строка - строка которую можно преобразовать в символьное выражение, используя eval
                                    # Вторая строка - строка, эквивалент первой, но записано в другой форме, для удобства мутаций
        self.get_eval_expr()        # <- Первая строка
                                    # 
        self.get_mut_expr()         # <- Вторая строка

    def get_eval_expr(self):
        self.eval_string = ""                           #
        self.get_eval_expr_rec(self.root)               # Рекурсивно получаем строку из дерева
        self.eval_string = self.eval_string.strip(' ')  # алгоритмы для обычной строки и для строки мутаций одинаковы

    def get_eval_expr_rec(self, root):
        """
        Рекурсивая функция получающая из дерева выражений
        человеческую строку выражения
        root: узел дерева 
        """
        # В строке для мутаций нет проверки на возведение в степень
        # так как, возведение в степень является бинарной операцией
        # ее можно также записать как и остальные в виде  ( * X * )
        # Но, для правильного перевода в числовое выражение удобнее использовать запись в виде функции X(*,*)
        if root is None:
            return
        # print(root)
        if(root.operation in self.functions):           # Если операция функция
            self.eval_string += " " + root.operation    # += | _func_(
            self.eval_string += " ("                    #
        
        if(root.operation in self.operators):           # Если операция
            self.eval_string += " ("                    # += | _( 


        if (root.operation == "apow"):                  # Если операция оказалась операцией возведения в степень
            self.eval_string += " apow ("               # += |_apow_(
            self.get_eval_expr_rec(root.left)           # += | _apow_(_left
            self.eval_string += " ,"                    # += | _apow_(_left_,
            self.get_eval_expr_rec(root.right)          # += | _apow_(_left_,_right
            self.eval_string += " )"                    # += | _apow_(_left_,_right_)
        else:
            self.get_eval_expr_rec(root.left)                           # += | _(_left
            if(not root.operation in self.functions):                   # Дело в том, что левый узел может представлять из себя тоже выражение,
                if(not root.operation in self.operators):               # но в крайнем случае представляет собой элемент, на всякий случай обособляю
                    self.eval_string += " " + root.operation            # их скобка (лишним не бывает)
                else:                                                   # 
                    self.eval_string += " " + root.operation            # += | _(_left[_(_symbol_)]
                                                                        # += | _(_left[_(_symbol_)]_oper
            self.get_eval_expr_rec(root.right)                          # += | _(_left[_(_symbol_)]_oper_right[_(_symbol_)] <- Нужна закрывающая скобка

        if(root.operation in self.operators):
            self.eval_string += " )" 

        if(root.operation in self.functions):
            self.eval_string += " )" 

    def get_mut_expr(self):
        # Для пояснений смотреть комментарии функции get_eval_expr
        self.mut_string = ""
        self.get_mut_expr_rec(self.root)             
        self.mut_string = self.mut_string.strip(' ')

    def get_mut_expr_rec(self, root):
        """
        Рекурсивая функция получающая из дерева выражений
        человеческую строку выражения
        root: узел дерева 
        """
        if root is None:
            return
        
        else:

            if(root.operation in self.functions):
                self.mut_string += " " + root.operation
                self.mut_string += " ("
            
            if(root.operation in self.operators):
                self.mut_string += " ("
            
            self.get_mut_expr_rec(root.left)
            if(not root.operation in self.functions):
                if(not root.operation in self.operators):
                    self.mut_string += " ( " + root.operation + " )"
                else:
                    self.mut_string += " " + root.operation
        
            self.get_mut_expr_rec(root.right)

            if(root.operation in self.operators):
                self.mut_string += " )" 

            if(root.operation in self.functions):
                self.mut_string += " )" 


    def mutations(self):
        
        for i in range(self.NOM):
            self.mutate_tree()
        
        self.get_eval_expr()
        self.get_mut_expr()

    def mutate_tree(self):
        self.mutate_tree_rec(self.root)
        self.tree_length = self.get_tree_length()

    def mutate_tree_rec(self, root):
        prob = np.random.rand()
        prob_add = np.random.rand()
        prob_add_un = np.random.rand()
        if(root == None):
            return
        
        # if(self.tree_length <= self.MLT):
        if(prob_add < 0.3):
            if(np.random.rand() > 0.3):
                x_str = "x" + str(np.random.randint(self.x_length))
                l = np.random.randint(len(self.operators))
                x_node = node(x_str)
                nod = node(self.operators[l], x_node, root)
                x_node.set_par(nod)

                if(root.par != None):
                    if(root.par.left == root):
                        root.par.left = nod
                        nod.set_par(root.par)
                    else:
                        root.par.right = nod
                        nod.set_par(root.par)
                else:
                    self.root = nod
                root.set_par(nod)

            else:
                const = 5*np.random.rand()
                l = np.random.randint(len(self.operators))
                c_node = node(str(const))
                nod = node(self.operators[l], c_node, root)
                c_node.par = nod

                if(root.par != None):
                    if(root.par.left == root):
                        root.par.left = nod
                        nod.set_par(root.par)
                    else:
                        root.par.right = nod
                        nod.set_par(root.par)
                else:
                    self.root = nod
                root.set_par(nod)

            self.mutate_tree_rec(root.right)
            self.mutate_tree_rec(root.left)
            return

        if(prob_add_un < 0.1):                             # Добавление унарной операции
            l = np.random.randint(len(self.functions))
            nod = node(self.functions[l], None, root)

            if(root.par != None):
                    if(root.par.left == root):
                        root.par.left = nod
                        nod.set_par(root.par)
                    else:
                        root.par.right = nod
                        nod.set_par(root.par)
            else:
                self.root = nod
            root.set_par(nod)

            self.mutate_tree_rec(root.right)
            self.mutate_tree_rec(root.left)
            return
        # else:
        #     prob -= 0.05
                
        if(root.operation in self.operators):
            if(prob < 0.1):                         # Замена бинарной операции
                l = np.random.randint(len(self.operators))
                root.operation = self.operators[l]
                self.mutate_tree_rec(root.right)
                self.mutate_tree_rec(root.left)
                return

        if(root.operation in self.functions):      
            if(prob < 0.1):                        # Удаление унарной операции
                if(root.par != None):
                    root.right.set_par(root.par)
                    if(root.par.right == root):

                        root.par.right = root.right
                    else:
                        root.par.left = root.right
                else:
                    root.right.set_par(root.par)
                    self.root = root.right

                root = root.right
                if(root.operation in self.functions):
                    self.mutate_tree_rec(root)
                    return
                self.mutate_tree_rec(root.right)
                self.mutate_tree_rec(root.left)
                return
                
            if(prob < 0.14 and prob > 0.07):                        # Замена унарной операции
                l = np.random.randint(len(self.functions))
                root.operation = self.functions[l]
                self.mutate_tree_rec(root.right)
                self.mutate_tree_rec(root.left)
                return

            # if(prob < 0.05):
            #     if(root.par != None):
            #         if(root.par.left == root):
            #             c_node = node(str(const))
            #             root.par.left = c_node
            #             c_node.set_par(root.par)
            #         else:
            #             c_node = node(str(const))
            #             root.par.right = c_node
            #             c_node.set_par(root.par)
            #         root.par = None

            #         return
                    
        
    def find_right_bracker_index(self, string, index):
        count = 0 

        for i in range(index + 1, len(string)):
            if(string[i] == "("):
                count += 1
            if(string[i] == ")"):
                count -= 1

            if(count == 0):
                return i
        
        return -1 

    def get_eval(self):
        eval_s = self.eval_string 
        return eval_s
    
    def get_mut(self):
        return self.mut_string

    def is_symbol_func(self,x):
        return (not (x in self.operators)) and (x in self.functions or (x != "(" and x != ")"))

    def get_tree_length(self):
        return self.get_tree_length_rec(self.root)

    def get_tree_length_rec(self, root):
        if(root == None):
            return 0

        if(root.right == None and root.left == None):
            return 0
        
        return 1 + max(self.get_tree_length_rec(root.right), self.get_tree_length_rec(root.left))

    def get_node_tree(self):
        return self.get_node_tree_rec(self.root)

    def get_node_tree_rec(self,node):
        if(node.left == None and node.right == None):
            return node
        
        prob = np.random.rand()
        if(prob < 0.1):
            return node

        if(node.left == None):
            return self.get_node_tree_rec(node.right)
        
        division = np.random.rand()
        if(division < 0.5):
            return  self.get_node_tree_rec(node.left)
        else: 
            return  self.get_node_tree_rec(node.right)
        
    def get_eval_by_tree(self, node):
        if node is None:
            return
        # print(root)
        if(node.operation in self.functions):
            print(" " + node.operation, end = "")
            print(" (", end = "")
        
        if(node.operation in self.operators):
            self.eval_string += " ("


        if (node.operation == "apow"):                  # Если операция оказалась операцией возведения в степень
            print(" apow (", end = "")              # += |_apow_(
            self.get_eval_by_tree(node.left)           # += | _apow_(_left
            print(" ,", end = "")                    # += | _apow_(_left_,
            self.get_eval_by_tree(node.right)          # += | _apow_(_left_,_right
            print(" )",end = "")                  # += | _apow_(_left_,_right_)
        else:
            self.get_eval_by_tree(node.left)                           # += | _(_left
            if(not node.operation in self.functions):                   # Дело в том, что левый узел может представлять из себя тоже выражение,
                if(not node.operation in self.operators):               # но в крайнем случае представляет собой элемент, на всякий случай обособляю
                    print(" " + node.operation + " ", end = "")   # их скобка (лишним не бывает)
                else:                                                   # 
                    print(" " + node.operation, end = "")            # += | _(_left[_(_symbol_)]
                                                                        # += | _(_left[_(_symbol_)]_oper
            self.get_eval_by_tree(node.right)                          # += | _(_left[_(_symbol_)]_oper_right[_(_symbol_)] <- Нужна закрывающая скобка

        if(node.operation in self.operators):
            print(" )", end = "") 

        if(node.operation in self.functions):
            print(" )", end = "")

    def get_evaluation_number(self, x_array):
        if(len(x_array) != self.x_length):
            return
        
        s = self.eval_string
        for i in range(self.x_length):
            x_str = "x" + str(i)
            s = s.replace(x_str, str(x_array[i]))
        
        return eval(s)

def alog(x):
    return log(Abs(x) + 0.00000001)

def apow(x,y):
    return (Pow(abs(x), y))


def get_score(expr: symb_exp,data_x, data_y):
        L = len(data_x[0])
        score = 0
        func = eval(expr.get_eval())
        print("OK")
        print(f"{expr.tree_length}")
        for i in range(L):
            x = []
            
            for j in range(3):      #
                x.append(data_x[j][i])     # Формирование вектора (точки)
            
            try:                                # Если число получается слишком большим
                if(score > 1e300):
                    1/0
                # score += (expr.get_evaluation_number(x) - data_y[i])**2/(L-1)
                score += func.subs({'x1' : x[0],'x1' : x[1],'x2': x[2]})
            except:                             # То представитель плохой сделаем копию лучшего
                return


        score = (score)**0.5
        expr.score = score
        print(score)

def test_function(x1,x2,x3):
    return 2*x1 + 10*np.log(x2) + np.exp(x3)

if __name__ == "__main__":
    x0, x1, x2, x3, x4, x5 = symbols('x:6')
    test0  = "sin ( alog ( 3.0 ) * 0.3 * x1 ) + alog ( apow ( 4 , 0.5 ) )"
    test1  = "(-1.0)* ( ( 4 ) ) apow ( 0.5 ) + 0.5"
    test2 = "sin ( ( 4 ) ) apow ( 0.5 ) apow ( 0.5 )"
    test3 = "alog ( x3 ) + apow ( x1 , 0.5 ) + sin ( x1 )"

    X = 10*np.random.rand(3,300)
    Y = test_function(X[0], X[1], X[2])

    exp1 = symb_exp(x_length = 3, number_of_mutations = 1)
    exp1.get_symb_str("1.0 + 1.0")
    
    r = [exp1]
    # exp1.mutations()
    for i in range(40):
        exp1.mutations()

    
    print(exp1.get_eval)


    # exp2 = symb_exp(x_length = 3, number_of_mutations = 1)
    # exp2.get_symb_str(exp1.get_eval())
    # print(exp2 == exp1)
    # print(exp2.get_eval() == exp1.get_eval())
    # print(exp2.get_eval())
    # print(exp1.get_eval())

    # print(eval(exp1.get_evaluation_number([1,1,1])))
    # print(eval(exp1.get_eval()).subs({'x0': 1,"x1" : 1 ,"x2" : 1}))
    # tree = exp1.get_node_tree()
    # exp1.eval_string = ""
    # exp1.get_eval_expr_rec(tree)
    # print(f"\n\n\n")
    # # exp1.get_eval_by_tree(exp1.root)
    
    # print(eval(exp1.get_eval()))
    # print(eval(exp1.get_mut()))

