PROJECT(RemoveStool4) 									
FIND_PACKAGE(ITK) 								
IF(ITK_FOUND) 									
	INCLUDE(${ITK_USE_FILE}) 					
ELSE(ITK_FOUND) 								
MESSAGE(FATAL_ERROR 							
"ITK not found. Please set ITK_DIR.") 			
ENDIF(ITK_FOUND) 								
 												
ADD_EXECUTABLE(RemoveStool4 RemoveStool4.cxx )  					
 
TARGET_LINK_LIBRARIES(RemoveStool4 ${ITK_LIBRARIES}) 	
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}) 		
INCLUDE_DIRECTORIES(C:/GitProjects/common) 	
