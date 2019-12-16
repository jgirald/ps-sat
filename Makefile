all : ps-sat

clean :
	rm ps-sat

%: %.cpp
	g++ $< -o $@
