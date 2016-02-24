__author__ = 'chengwei'
import os
import time
import argparse

class hw3():
    def __init__(self,filename):
        (self.inputsize,self.inputline,self.file)=self.read_file(filename)
        self.filename = filename

    def divide_conquer(self,list,left,right):
        if left == right:
            return (left,left, list[left])
        mid = (left+right )/2
        (l_i,l_j,left_sum) = self.divide_conquer(list,left,mid)
        (r_i,r_j,right_sum) = self.divide_conquer(list,mid+1,right)
        (m_i,m_j,m_sum)=self.mid(list,left,mid,right)
        total=max(left_sum,m_sum,right_sum)
        if total == left_sum:
            return (l_i,l_j,left_sum)
        elif total == right_sum:
            return (r_i,r_j,right_sum)
        else:
            return (m_i,m_j,m_sum)
    def mid(self,list,left,mid,right):
        i=mid
        sum=0.0
        lef_r =-float('inf')
        for k in range(mid,left-1,-1):
            sum+=list[k]
            if lef_r < sum:
                i=k
                lef_r=sum
        j=mid+1;rig_r = -float('inf');sum=0
        for k in range(mid+1,right+1):
            sum+=list[k]
            if sum>rig_r:
                j=k;rig_r=sum
        total = rig_r+lef_r
        return (i,j,total)

    def dynamic(self,list):
        if len(list) == 0 or len(list) ==1:
            return list;
        memory=[0]*len(list)
        memory[0]=max(list[0],0)
        for i in range(1,len(list)):
            memory[i]=max(memory[i-1]+list[i],0)
        #back tracking
        get_max=lambda s: max((x, i) for i, x in enumerate(s))[1]
        j=get_max(memory)
        i=j
        while memory[i-1] >0 and i >0 :
            i-=1
        return (memory[j],i+1,j+1)

    def read_file(self,filename):
        a= os.path.dirname(os.path.realpath(__file__))
        a=a+'/data/'+filename
        size=0
        line_no=0
        input=[]
        with open(a) as file:
            lines = file.readlines();
            first= [int(e.strip()) for e in lines[0].split(',')]
            size=first[0]
            line_no=first[1]
            for each_input in lines[1:]:
                input.append([float(e.strip()) for   e in each_input.split(',')])
        return (size,line_no,input)
    def run_output(self):
        path=os.path.dirname(os.path.realpath(__file__))+'/output/'
        output=open(path+'cli402_output_dc_'+str(self.inputsize)+'.txt','w')
        for i in range(self.inputline):
            begin=time.time()
            result=self.divide_conquer(self.file[i],0,self.inputsize-1)
            dura=(time.time()-begin)*1000
            (i,j,sum)=result;
            sum=float("%.2f" %sum)
            result=(sum,i+1,j+1)
            string_result='{}'.format((result,dura)).replace('(','').replace(')','')+'\n'
            output.write(string_result)
        output.close()

        output_dp=open(path+'cli402_output_dp_'+str(self.inputsize)+'.txt','w')
        for j in range(self.inputline):
            begin=time.time()
            result=self.dynamic(self.file[j])
            dura=(time.time()-begin)*1000
            (sum,i,j)=result
            sum=float("%.2f" %sum)
            result=(sum,i,j)
            string_result='{}'.format((result,dura)).replace('(','').replace(')','')+'\n'
            output_dp.write(string_result)
        output_dp.close()

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="input data name")
    parser.add_argument('name',help='just input the name of the input file,no path attached')
    args=parser.parse_args()
    a=hw3(args.name)
    a.run_output()