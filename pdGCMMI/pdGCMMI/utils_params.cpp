/*
version 1.1
 -  key���� �������� ������ �߻��Ǵ� ���� ����

version 1.0
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "utils_params.h"

/**
* @fn int getProfileString(FILE* fp, char* section, char* key, char* pBuf, int bufLen)
* @brief ���Ͽ��� section ī������ key���� �б�
* @param fp ����
* @param section section ����
* @param key key ����
* @param pBuf key���� ����� ���ڿ� ������ bufLen�� ũ�⸦ ���´�.
* @param bufLen key�� ���ڿ��� �ִ� ����
* @return -1 �Ķ���� ����
* @return -2 ���� �ڵ� ����
* @return -3 section �˻� ����
* @return -4 key �˻� ����
* @return 0 ����
 */

int getProfileString(FILE* fp, const char* section, const char* key, char* pBuf, int bufLen) {

	char bufKey[MAX_LENGTH];
	char bufValue[MAX_LENGTH];
	int checkP;

	// section, key �Ķ���� ���� üũ
	if ((strlen(key) < 1) || (strlen(key) < 1)) return -1;

	memset(pBuf, 0, sizeof(char) * bufLen);

	//���� �ڵ� ����
	if (fp == NULL) return -2;

	// fp�� ��ġ�� ������ ���� ������ �̵�
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

	// Key���� ������
	if (checkP == 0) {
		pBuf[0] = '\0';
		return -4;
	}

	return 0;
}


/**
* @fn int getNextKey(FILE* fp, char* retKey, char* retValue)
* @brief ���Ͽ��� key���� ���� ������ �����Ͽ� retKey�� retValue�� ����
* @param fp ����
* @param retKey key�� ī�װ� ����
* @param retValue key�� ����
* @return -1 �Ķ���� ����
* @return -2 ���� �ڵ� ����
* @return -3 section ����
* @return 0 ����
 */
int getNextKey(FILE* fp, char* retKey, char* retValue) {

	char str[MAX_LENGTH];
	char buf[MAX_LENGTH];
	char* p = NULL;
	int checkP;
	int size;

	// �Ķ���� ���� üũ
	if (retKey == NULL || retValue == NULL) return -1;

	//���� �ڵ� ����
	if (fp == NULL) return -2;

	// serach Key
	checkP = 0;

	// feof() ������ ���̸� 0�� �ƴ� �� ����
	while (!feof(fp)) {
		fgets(buf, MAX_LENGTH, fp);
		trimString(buf);
		if (buf[0] == 0 || isComment(buf[0])) continue;
		if (buf[0] == '[') return -3;

		break;
	}

	// token�� �˻�
	p = buf;
	size = 0;
	while ((*p != '\n') && (*p != '=') && (*p)) {
		size++;
		p++;
	}

	// p�� '=' �̸� value�� �ִ°���
	if (*p == '=') checkP = 1;
	*p = '\0';

	trimString(buf);

	// version 1.1 �������
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
* @brief ���Ͽ��� section�� int���� key ���� �д´�.
* @param fp ����
* @param section �˻��� section �̸�
* @param key �˻��� key �̸�
* @param ret key���� ����� �ּ�
* @return -1 �Ķ���� ����
* @return -2 ���� �ڵ� ����
* @return -3 section �˻� ����
* @return -4 key �˻� ����
* @return 0 ����
 */
int getIniInt(FILE* fp, const char* section, const char* key, int* ret) {
	char buf[MAX_LENGTH];
	int n = 0;

	// �Ķ���� ����
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
* @brief ���Ͽ��� section�� float���� key ���� �д´�.
* @param fp ����
* @param section �˻��� section �̸�
* @param key �˻��� key �̸�
* @param ret key���� ����� �ּ�
* @return -1 �Ķ���� ����
* @return -2 ���� �ڵ� ����
* @return -3 section �˻� ����
* @return -4 key �˻� ����
* @return 0 ����
 */
int getIniFloat(FILE* fp, const char* section, const char* key, float* ret) {
	char buf[MAX_LENGTH];
	int n = 0;

	// �Ķ���� ����
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
* @brief fp�� ���Ͽ��� section ī�װ��� ���������� �̵���Ų��.
* @param fp ����
* @param section �˻��� section �̸�
* @return -1 �Ķ���� ����
* @return -2 ���� �ڵ� ����
* @return -3 section �˻� ����
* @return 0 ����
 */
int seekSection(FILE* fp, const char* section) {

	char str[MAX_LENGTH];
	char buf[MAX_LENGTH];

	int checkP;

	//�Ķ���� ���� üũ
	if (strlen(section) < 1) return -1;

	// ���� �ڵ� ���� üũ
	if (fp == NULL) return -2;

	sprintf_s(str, "[%s]", section);

	// search key
	checkP = 0;
	// feof() ������ ���̸� 0�� �ƴ� �� ����
	while (!feof(fp)) {
		// ���Ͽ��� ������ �о��
		fgets(buf, MAX_LENGTH, fp);

		// ;�� �ּ���
		if (buf[0] == 0 || isComment(buf[0])) continue;

		trimString(buf);

		if (strncmp(buf, str, strlen(str)) == 0) {
			checkP = 1;
			break;
		}
	}

	// �˻��� Section�� ���� ���
	if (checkP == 0)
		return -3;

	return 0;
}

/**
* @fn char* trimString(char* str)
* @brief �Էµ� ���ڿ����� ��, �� ������ ����
* @details  ���ڿ��� �߰��� �߻��� ��ĭ�� ���� ��ȯ������ ������ ���ŵ� ���ڿ��� ������, �Էµ� str���� �����
* @param str ���ڿ�
* @return str ������ ���ŵ� ���ڿ�
 */

char* trimString(char* str) {

	char* f = str, * e = NULL, * c = str;

	// ���� ���� ����
	e = str + (strlen(str)) - 1;

	while (isBlank(*e) && f <= e) e--;
	*(e + 1) = '\0';

	//���� ���� ����
	while (isBlank(*f) && f <= e) f++;

	//������� �κ� ����
	if (str != f) {
		while (f <= e) *(c++) = *(f++);

		*c = '\0';
	}

	return str;
}


/**
* @fn bool isBlank(char ch)
* @brief �Էµ� ���ڰ� �������� �Ǵ��ϴ� �Լ�
* @param ch ����
* @return True ���� ������ ���
* @return False ������ �ƴ� ���
*/
bool isBlank(char ch) {
	// 0x20: Space
	return ((ch == 0x20) || (ch == '\t') || (ch == '\r') || (ch == '\n'));
}

/**
* @fn isComment(char ch)
* @brief �Էµ� ���ڰ� �ּ��� �������� �Ǵ�
* @param ch ����
* return True �ּ� ����
* return False �ּ� ���� �ƴ�
*/
bool isComment(char ch) {
	return ((ch == ';'));
}