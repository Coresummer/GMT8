MCL_DIR?=../mcl
CXX = clang++-12
PROGRAM = CP8

OBJS   =  main.o create.o scalar.o mpn.o fp.o fp2.o field_test.o time.o fp4.o fp8.o efp.o efp2.o efp4.o efp8.o test_efp.o miller.o final_exp.o test_pairing.o mcl.o 
HEADER = define.h create.h  scalar.h mpn.h fp.h fp2.h fp2.h field_test.h time.h fp4.h fp8.h efp.h efp2.h efp4.h efp8.h test_efp.h miller.h final_exp.h test_pairing.h 

$(PROGRAM): $(OBJS)
	$(CXX) -fPIC -g -pg -o $(PROGRAM) $(OBJS) -Ofast -march=native -lgmp -lstdc++
	#-fsanitize=address

#サフィックスルールの適応対象の拡張子の定義
.SUFFIXES: .c .o

CFLAGS+=-Ofast -fPIC -g -DNDEBUG -I $(MCL_DIR)/include -I $(MCL_DIR)/src
#c言語のソースファイルから拡張子が.oのファイルを作成する

.c.o:
	$(CC) -c $< $(CFLAGS) -MMD -MP -MF $(@:.o=.d)
.cpp.o:
	$(CXX) -c $< $(CFLAGS) -MMD -MP -MF $(@:.o=.d)

#ヘッダファイルの依存関係(ヘッダファイルが変わったらすべてコンパイルし直す)
$(OBJS): $(HEADER)

DEPS=$(OBJS:.o=.d)
-include $(DEPS)

#不要なファイル削除用(コマンド:make clean)
.PHONY: clean
clean:
	rm -f $(PROGRAM) *.o *.d

# don't remove these files automatically
.SECONDARY: $(OBJS)
