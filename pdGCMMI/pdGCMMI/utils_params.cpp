/*
version 1.1
 -  key값이 존재하지 않을때 발생되는 오류 수정

version 1.0
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils_params.h"

/**
* @fn int getProfileString(FILE* fp, char* section, char* key, char* pBuf, int bufLen)
* @brief 파일에서 section 카데고리의 key값을 읽기
* @param fp 파일
* @param section section 정보
* @param key key 정보
* @param pBuf key값이 저장될 문자열 변수로 bufLen의 크기를 갖는다.
* @param bufLen key값 문자열의 최대 길이
* @return -1 파라미터 에러
* @return -2 파일 핸들 에러
* @return -3 section 검색 실패
* @return -4 key 검색 실패
* @return 0 성공
 */

int getProfileString(FILE* fp, const char* section, const char* key, char* pBuf, int bufLen) {

	char bufKey[MAX_LENGTH];
	char bufValue[MAX_LENGTH];
	int checkP;

	// section, key 파라미터 오류 체크
	if ((strlen(key) < 1) || (strlen(key) < 1)) return -1;

	memset(pBuf, 0, sizeof(char) * bufLen);

	//파일 핸들 에러
	if (fp == NULL) return -2;

	// fp의 위치를 파일의 가장 앞으로 이동
	fseek(fp, 0L, SEEK_SET);

	if (seekSection(fp, section) != 0) return -3;

	// search Key
	getNextKey(fp, bufKey, bufValue);

	checkP = 0;
	while (bufKey[0]) {
		if (strcmp(bufKey, key) == 0) {
			strncpy_s(pBuf, bufLen, bufValue, strlen(bufValue));
			pBuf[strlen(bufValue) + 1] = '\0';
			checkP = 1;

			break;
		}
		if (getNextKey(fp, bufKey, bufValue))
			break;
	}

	// Key값이 없으면
	if (checkP == 0) {
		pBuf[0] = '\0';
		return -4;
	}

	return 0;
}


/**
* @fn int getNextKey(FILE* fp, char* retKey, char* retValue)
* @brief 파일에서 key값에 대한 정보를 추출하여 retKey와 retValue에 저장
* @param fp 파일
* @param retKey key값 카테고리 저장
* @param retValue key값 저장
* @return -1 파라미터 에러
* @return -2 파일 핸들 에러
* @return -3 section 종료
* @return 0 성공
 */
int getNextKey(FILE* fp, char* retKey, char* retValue) {

	char str[MAX_LENGTH];
	char buf[MAX_LENGTH];
	char* p = NULL;
	int checkP;
	int size;

	// 파라미터 오류 체크
	if (retKey == NULL || retValue == NULL) return -1;

	//파일 핸들 에러
	if (fp == NULL) return -2;

	// serach Key
	checkP = 0;

	// feof() 파일의 끝이면 0이 아닌 수 리턴
	while (!feof(fp)) {
		fgets(buf, MAX_LENGTH, fp);
		trimString(buf);
		if (buf[0] == 0 || isComment(buf[0])) continue;
		if (buf[0] == '[') return -3;

		break;
	}

	// token을 검색
	p = buf;
	size = 0;
	while ((*p != '\n') && (*p != '=') && (*p)) {
		size++;
		p++;
	}

	// p가 '=' 이면 value가 있는것임
	if (*p == '=') checkP = 1;
	*p = '\0';

	trimString(buf);

	// version 1.1 변경사항
	if (strlen(buf) < MAX_LENGTH) {
		strncpy_s(retKey, MAX_LENGTH, buf, strlen(buf));
		retKey[strlen(buf) + 1] = '\0';
	}
	else {
		strncpy_s(retKey, MAX_LENGTH, buf, MAX_LENGTH - 2);
		retKey[MAX_LENGTH - 1] = '\0';
		return -4;
	}

	// value
	if (checkP == 1) {
		strcpy_s(str, p + 1);

		p = str;
		while ((*p != '\n') && (*p)) p++;
		*p = '\0';

		trimString(str);

		strncpy_s(retValue, MAX_LENGTH, str, strlen(str));
		retValue[strlen(str) + 1] = '\0';
	}

	return 0;
}

/**
* @fn int getIniInt(FILE* fp, char* section, char* key, int* ret)
* @brief 파일에서 section의 int형의 key 값을 읽는다.
* @param fp 파일
* @param section 검색할 section 이름
* @param key 검색할 key 이름
* @param ret key값이 저장될 주소
* @return -1 파라미터 에러
* @return -2 파일 핸들 에러
* @return -3 section 검색 실패
* @return -4 key 검색 실패
* @return 0 성공
 */
int getIniInt(FILE* fp, const char* section, const char* key, int* ret) {
	char buf[MAX_LENGTH];
	int n = 0;

	// 파라미터 에러
	if (ret == 0) return -1;

	n = getProfileString(fp, section, key, buf, MAX_LENGTH);

	if (n == 0)
	{
		(*ret) = atoi(buf);
		return 0;
	}

	return n;
}


/**
* @fn int getIniFloat(FILE* fp, char* section, char* key, float* ret)
* @brief 파일에서 section의 float형의 key 값을 읽는다.
* @param fp 파일
* @param section 검색할 section 이름
* @param key 검색할 key 이름
* @param ret key값이 저장될 주소
* @return -1 파라미터 에러
* @return -2 파일 핸들 에러
* @return -3 section 검색 실패
* @return -4 key 검색 실패
* @return 0 성공
 */
int getIniFloat(FILE* fp, const char* section, const char* key, float* ret) {
	char buf[MAX_LENGTH];
	int n = 0;

	// 파라미터 에러
	if (ret == 0)
		return -1;

	n = getProfileString(fp, section, key, buf, MAX_LENGTH);

	if (n == 0)
	{
		(*ret) = atof(buf);
		return 0;
	}

	return n;
}


/**
* @fn int seekSection(FILE* fp, char* section)
* @brief fp를 파일에서 section 카테고리의 시작점으로 이동시킨다.
* @param fp 파일
* @param section 검색할 section 이름
* @return -1 파라미터 에러
* @return -2 파일 핸들 에러
* @return -3 section 검색 실패
* @return 0 성공
 */
int seekSection(FILE* fp, const char* section) {

	char str[MAX_LENGTH];
	char buf[MAX_LENGTH];

	int checkP;

	//파라미터 에러 체크
	if (strlen(section) < 1) return -1;

	// 파일 핸들 에러 체크
	if (fp == NULL) return -2;

	sprintf_s(str, "[%s]", section);

	// search key
	checkP = 0;
	// feof() 파일의 끝이면 0이 아닌 수 리턴
	while (!feof(fp)) {
		// 파일에서 한줄을 읽어옴
		fgets(buf, MAX_LENGTH, fp);

		// ;은 주석임
		if (buf[0] == 0 || isComment(buf[0])) continue;

		trimString(buf);

		if (strncmp(buf, str, strlen(str)) == 0) {
			checkP = 1;
			break;
		}
	}

	// 검색된 Section이 없을 경우
	if (checkP == 0)
		return -3;

	return 0;
}

/**
* @fn char* trimString(char* str)
* @brief 입력된 문자열에서 앞, 뒤 공백을 제거
* @details  문자열의 중간에 발생된 빈칸은 무시 반환값으로 공백이 제거된 문자열을 보내고, 입력된 str에도 적용됨
* @param str 문자열
* @return str 공백이 제거된 문자열
 */

char* trimString(char* str) {

	char* f = str, * e = NULL, * c = str;

	// 뒤쪽 공백 제거
	e = str + (strlen(str)) - 1;

	while (isBlank(*e) && f <= e) e--;
	*(e + 1) = '\0';

	//앞쪽 공백 제거
	while (isBlank(*f) && f <= e) f++;

	//공백없는 부분 복사
	if (str != f) {
		while (f <= e) *(c++) = *(f++);

		*c = '\0';
	}

	return str;
}


/**
* @fn bool isBlank(char ch)
* @brief 입력된 문자가 공백인지 판단하는 함수
* @param ch 문자
* @return True 공백 문자일 경우
* @return False 공백이 아닐 경우
*/
bool isBlank(char ch) {
	// 0x20: Space
	return ((ch == 0x20) || (ch == '\t') || (ch == '\r') || (ch == '\n'));
}

/**
* @fn isComment(char ch)
* @brief 입력된 문자가 주석의 시작인지 판단
* @param ch 문자
* return True 주석 문자
* return False 주석 문자 아님
*/
bool isComment(char ch) {
	return ((ch == ';'));
}